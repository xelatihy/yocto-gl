//
// Yocto/Draw: Utilities for real-time reandering of a scene.
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

#ifndef _YOCTO_DRAW_
#define _YOCTO_DRAW_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <string>
#include <vector>

#include "yocto_opengl.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

}  // namespace yocto

namespace yocto {

// Opengl caemra
struct gui_camera {
  frame3f frame    = identity3x4f;
  float   lens     = 0.050;
  float   aspect   = 1.000;
  float   film     = 0.036;
  float   near     = 0.001;
  float   far      = 10000;
  float   aperture = 0;
  float   focus    = 0;
};

// Opengl texture
struct gui_texture {
  // shape properties
  ogl_texture* texture = new ogl_texture{};

  // Disable copy construction
  gui_texture()                   = default;
  gui_texture(const gui_texture&) = delete;
  gui_texture& operator=(const gui_texture&) = delete;

  // Cleanup
  ~gui_texture();
};

// Opengl material
struct gui_material {
  // material
  vec3f        emission      = {0, 0, 0};
  vec3f        color         = {0, 0, 0};
  float        metallic      = 0;
  float        roughness     = 0;
  float        specular      = 0;
  float        opacity       = 1;
  gui_texture* emission_tex  = nullptr;
  gui_texture* color_tex     = nullptr;
  gui_texture* metallic_tex  = nullptr;
  gui_texture* roughness_tex = nullptr;
  gui_texture* specular_tex  = nullptr;
  gui_texture* opacity_tex   = nullptr;
  gui_texture* normal_tex    = nullptr;
};

struct gui_shape {
  // shape properties
  ogl_shape* shape = new ogl_shape{};

  // Disable copy construction
  gui_shape()                 = default;
  gui_shape(const gui_shape&) = delete;
  gui_shape& operator=(const gui_shape&) = delete;

  // Cleanup
  ~gui_shape();
};

// Shading type
enum struct gui_shading_type { constant = 0, shaded };

// Opengl instance
struct gui_instance {
  // instance properties
  frame3f          frame       = identity3x4f;
  gui_shape*       shape       = nullptr;
  gui_material*    material    = nullptr;
  bool             hidden      = false;
  bool             highlighted = false;
  gui_shading_type shading     = gui_shading_type::shaded;
};

// Opengl environment
struct gui_environment {
  // environment properties
  frame3f      frame        = identity3x4f;
  vec3f        emission     = {1, 1, 1};
  gui_texture* emission_tex = nullptr;

  // drawing data
  gui_shape*   shape   = new gui_shape{};
  ogl_cubemap* cubemap = new ogl_cubemap{};

  // envlight precomputed data
  ogl_cubemap* envlight_diffuse  = new ogl_cubemap{};
  ogl_cubemap* envlight_specular = new ogl_cubemap{};
  ogl_texture* envlight_brdflut  = new ogl_texture{};

  // Disable copy construction
  gui_environment()                       = default;
  gui_environment(const gui_environment&) = delete;
  gui_environment& operator=(const gui_environment&) = delete;

  // Cleanup
  ~gui_environment();
};

// Opengl scene
struct gui_scene {
  // scene objects
  vector<gui_camera*>      cameras      = {};
  vector<gui_instance*>    instances    = {};
  vector<gui_shape*>       shapes       = {};
  vector<gui_material*>    materials    = {};
  vector<gui_texture*>     textures     = {};
  vector<gui_environment*> environments = {};

  // programs
  ogl_program* environment_program = new ogl_program{};
  ogl_program* camlight_program    = new ogl_program{};
  ogl_program* envlight_program    = new ogl_program{};

  // disable copy construction
  gui_scene()                 = default;
  gui_scene(const gui_scene&) = delete;
  gui_scene& operator=(const gui_scene&) = delete;

  // cleanup
  ~gui_scene();
};

// Shading type
enum struct gui_lighting_type {
  envlight,
  camlight,
};

// Shading name
const auto gui_lighting_names = vector<string>{"envlight", "camlight"};

// Draw options
struct gui_scene_params {
  int               resolution       = 1280;
  bool              wireframe        = false;
  gui_lighting_type lighting         = gui_lighting_type::camlight;
  float             exposure         = 0;
  float             gamma            = 2.2f;
  bool              faceted          = false;
  bool              double_sided     = true;
  bool              non_rigid_frames = true;
  float             near             = 0.01f;
  float             far              = 10000.0f;
  vec4f             background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void init_scene(gui_scene* scene);
bool is_initialized(const gui_scene* scene);

// Initialize data for environment lighting
void init_environments(gui_scene* scene, bool precompute_envlight = true);

// Clear an OpenGL scene
void clear_scene(gui_scene* scene);

// old interface
[[deprecated]] void init_scene(gui_scene* scene, gui_texture* environment_tex,
    const vec3f& environment_emission);

// add scene elements
gui_camera*      add_camera(gui_scene* scene);
gui_texture*     add_texture(gui_scene* scene);
gui_material*    add_material(gui_scene* scene);
gui_shape*       add_shape(gui_scene* scene);
gui_instance*    add_instance(gui_scene* scene);
gui_environment* add_environment(gui_scene* scene);

// camera properties
void set_frame(gui_camera* camera, const frame3f& frame);
void set_lens(gui_camera* camera, float lens, float aspect, float film);
void set_nearfar(gui_camera* camera, float near, float far);

// check if initialized
bool is_initialized(const gui_texture* texture);
// clear texture
void clear_texture(gui_texture* texture);

// set texture
void set_texture(gui_texture* texture, const image<vec4b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui_texture* texture, const image<vec4f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(gui_texture* texture, const image<vec3b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui_texture* texture, const image<vec3f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(gui_texture* texture, const image<byte>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui_texture* texture, const image<float>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);

// material properties
void set_emission(gui_material* material, const vec3f& emission,
    gui_texture* emission_tex = nullptr);
void set_color(gui_material* material, const vec3f& color,
    gui_texture* color_tex = nullptr);
void set_metallic(gui_material* material, float metallic,
    gui_texture* metallic_tex = nullptr);
void set_roughness(gui_material* material, float roughness,
    gui_texture* roughness_tex = nullptr);
void set_specular(gui_material* material, float specular,
    gui_texture* specular_tex = nullptr);
void set_opacity(
    gui_material* material, float opacity, gui_texture* opacity_tex = nullptr);
void set_normalmap(gui_material* material, gui_texture* normal_tex);

// cheeck if initialized
bool is_initialized(const gui_shape* shape);
// clear
void clear_shape(gui_shape* shape);

// shape properties
void set_points(gui_shape* shape, const vector<int>& points);
void set_lines(gui_shape* shape, const vector<vec2i>& lines);
void set_triangles(gui_shape* shape, const vector<vec3i>& triangles);
void set_quads(gui_shape* shape, const vector<vec4i>& quads);
void set_positions(gui_shape* shape, const vector<vec3f>& positions);
void set_normals(gui_shape* shape, const vector<vec3f>& normals);
void set_texcoords(gui_shape* shape, const vector<vec2f>& texcoords);
void set_colors(gui_shape* shape, const vector<vec4f>& colors);
void set_tangents(gui_shape* shape, const vector<vec4f>& tangents);
void set_instance_from(gui_shape* shape, const vector<vec3f>& froms);
void set_instance_to(gui_shape* shape, const vector<vec3f>& tos);

// get shaoe properties
ogl_arraybuffer* get_positions(gui_shape* shape);
ogl_arraybuffer* get_normals(gui_shape* shape);
ogl_arraybuffer* get_texcoords(gui_shape* shape);
ogl_arraybuffer* get_colors(gui_shape* shape);
ogl_arraybuffer* get_tangents(gui_shape* shape);

// instance properties
void set_frame(gui_instance* instance, const frame3f& frame);
void set_shape(gui_instance* instance, gui_shape* shape);
void set_material(gui_instance* instance, gui_material* material);
void set_hidden(gui_instance* instance, bool hidden);
void set_highlighted(gui_instance* instance, bool highlighted);

// check if initialized
bool is_initialized(const gui_environment* environment);
// clear environment
void clear_environment(gui_environment* environment);

// environment properties
void set_frame(gui_environment* environment, const frame3f& frame);
void set_emission(gui_environment* environment, const vec3f& emission,
    gui_texture* emission_tex = nullptr);

// shortcuts
gui_camera*      add_camera(gui_scene* scene, const frame3f& frame, float lens,
         float aspect, float film = 0.036, float near = 0.001, float far = 10000);
gui_material*    add_material(gui_scene* scene, const vec3f& emission,
       const vec3f& color, float specular, float metallic, float roughness,
       gui_texture* emission_tex = nullptr, gui_texture* color_tex = nullptr,
       gui_texture* specular_tex = nullptr, gui_texture* metallic_tex = nullptr,
       gui_texture* roughness_tex = nullptr, gui_texture* normalmap_tex = nullptr);
gui_shape*       add_shape(gui_scene* scene, const vector<int>& points,
          const vector<vec2i>& lines, const vector<vec3i>& triangles,
          const vector<vec4i>& quads, const vector<vec3f>& positions,
          const vector<vec3f>& normals, const vector<vec2f>& texcoords,
          const vector<vec4f>& colors, bool edges = false);
gui_instance*    add_instance(gui_scene* scene, const frame3f& frame,
       gui_shape* shape, gui_material* material, bool hidden = false,
       bool highlighted = false);
gui_environment* add_environment(gui_scene* scene, const frame3f& frame,
    const vec3f& emission, gui_texture* emission_tex = nullptr);

// internal drawing functions
struct gui_scene_view {
  frame3f          camera_frame      = {};
  mat4f            view_matrix       = {};
  mat4f            projection_matrix = {};
  gui_scene_params params            = {};
};

void set_scene_view_uniforms(ogl_program* program, const gui_scene_view& view);
void set_instance_uniforms(ogl_program* program, const frame3f& frame,
    const gui_shape* shape, const gui_material* material, int shading_type = 0,
    bool double_sided = true, bool non_rigid_frames = false);
void set_eyelight_uniforms(ogl_program* program, const gui_scene_view& view);
void set_ibl_uniforms(ogl_program* program, const gui_scene* scene);

void draw_instances(gui_scene* scene, const gui_scene_view& view);
void draw_environment(gui_scene* scene, const gui_scene_view& view);

void draw_scene(gui_scene* scene, gui_camera* camera, const vec4i& viewport,
    const gui_scene_params& params);

// read-only access to defualt shader code
const char* draw_instances_vertex_code();
const char* draw_instances_eyelight_fragment_code();
const char* draw_instances_ibl_fragment_code();
const char* draw_enivronment_fragment_code();

}  // namespace yocto

#endif
