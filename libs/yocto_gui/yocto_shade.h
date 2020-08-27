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
struct shade_camera {
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
struct shade_texture {
  // shape properties
  ogl_texture* texture = new ogl_texture{};

  // Disable copy construction
  shade_texture()                     = default;
  shade_texture(const shade_texture&) = delete;
  shade_texture& operator=(const shade_texture&) = delete;

  // Cleanup
  ~shade_texture();
};

// Opengl material
struct shade_material {
  // material
  vec3f          emission      = {0, 0, 0};
  vec3f          color         = {0, 0, 0};
  float          metallic      = 0;
  float          roughness     = 0;
  float          specular      = 0;
  float          opacity       = 1;
  shade_texture* emission_tex  = nullptr;
  shade_texture* color_tex     = nullptr;
  shade_texture* metallic_tex  = nullptr;
  shade_texture* roughness_tex = nullptr;
  shade_texture* specular_tex  = nullptr;
  shade_texture* opacity_tex   = nullptr;
  shade_texture* normal_tex    = nullptr;
  bool           unlit         = false;
};

struct shade_shape {
  // shape properties
  ogl_shape* shape = new ogl_shape{};

  // Disable copy construction
  shade_shape()                   = default;
  shade_shape(const shade_shape&) = delete;
  shade_shape& operator=(const shade_shape&) = delete;

  // Cleanup
  ~shade_shape();
};

// Shading type
enum struct shade_shading_type { constant = 0, shaded };

// Opengl instance
struct shade_instance {
  // instance properties
  frame3f         frame       = identity3x4f;
  shade_shape*    shape       = nullptr;
  shade_material* material    = nullptr;
  bool            hidden      = false;
  bool            highlighted = false;
};

// Opengl environment
struct shade_environment {
  // environment properties
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {1, 1, 1};
  shade_texture* emission_tex = nullptr;

  // drawing data
  shade_shape* shape   = new shade_shape{};
  ogl_cubemap* cubemap = new ogl_cubemap{};

  // envlight precomputed data
  ogl_cubemap* envlight_diffuse  = new ogl_cubemap{};
  ogl_cubemap* envlight_specular = new ogl_cubemap{};
  ogl_texture* envlight_brdflut  = new ogl_texture{};

  // Disable copy construction
  shade_environment()                         = default;
  shade_environment(const shade_environment&) = delete;
  shade_environment& operator=(const shade_environment&) = delete;

  // Cleanup
  ~shade_environment();
};

// Opengl scene
struct shade_scene {
  // scene objects
  vector<shade_camera*>      cameras      = {};
  vector<shade_instance*>    instances    = {};
  vector<shade_shape*>       shapes       = {};
  vector<shade_material*>    materials    = {};
  vector<shade_texture*>     textures     = {};
  vector<shade_environment*> environments = {};

  // programs
  ogl_program* environment_program = new ogl_program{};
  ogl_program* camlight_program    = new ogl_program{};
  ogl_program* envlight_program    = new ogl_program{};

  // disable copy construction
  shade_scene()                   = default;
  shade_scene(const shade_scene&) = delete;
  shade_scene& operator=(const shade_scene&) = delete;

  // cleanup
  ~shade_scene();
};

// Shading type
enum struct shade_lighting_type {
  envlight,
  camlight,
};

// Shading name
const auto shade_lighting_names = vector<string>{"envlight", "camlight"};

// Draw options
struct shade_params {
  int                 resolution       = 1280;
  bool                wireframe        = false;
  shade_lighting_type lighting         = shade_lighting_type::camlight;
  float               exposure         = 0;
  float               gamma            = 2.2f;
  bool                faceted          = false;
  bool                double_sided     = true;
  bool                non_rigid_frames = true;
  float               near             = 0.01f;
  float               far              = 10000.0f;
  vec4f               background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void init_scene(shade_scene* scene);
bool is_initialized(const shade_scene* scene);

// Initialize data for environment lighting
void init_environments(shade_scene* scene, bool precompute_envlight = true);

// Check if we have an envlight
bool has_envlight(const shade_scene* scene);

// Clear an OpenGL scene
void clear_scene(shade_scene* scene);

// add scene elements
shade_camera*      add_camera(shade_scene* scene);
shade_texture*     add_texture(shade_scene* scene);
shade_material*    add_material(shade_scene* scene);
shade_shape*       add_shape(shade_scene* scene);
shade_instance*    add_instance(shade_scene* scene);
shade_environment* add_environment(shade_scene* scene);

// camera properties
void set_frame(shade_camera* camera, const frame3f& frame);
void set_lens(shade_camera* camera, float lens, float aspect, float film);
void set_nearfar(shade_camera* camera, float near, float far);

// check if initialized
bool is_initialized(const shade_texture* texture);
// clear texture
void clear_texture(shade_texture* texture);

// set texture
void set_texture(shade_texture* texture, const image<vec4b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(shade_texture* texture, const image<vec4f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(shade_texture* texture, const image<vec3b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(shade_texture* texture, const image<vec3f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(shade_texture* texture, const image<byte>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(shade_texture* texture, const image<float>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);

// material properties
void set_emission(shade_material* material, const vec3f& emission,
    shade_texture* emission_tex = nullptr);
void set_color(shade_material* material, const vec3f& color,
    shade_texture* color_tex = nullptr);
void set_metallic(shade_material* material, float metallic,
    shade_texture* metallic_tex = nullptr);
void set_roughness(shade_material* material, float roughness,
    shade_texture* roughness_tex = nullptr);
void set_specular(shade_material* material, float specular,
    shade_texture* specular_tex = nullptr);
void set_opacity(shade_material* material, float opacity,
    shade_texture* opacity_tex = nullptr);
void set_normalmap(shade_material* material, shade_texture* normal_tex);
void set_unlit(shade_material* material, bool unlit);

// cheeck if initialized
bool is_initialized(const shade_shape* shape);
// clear
void clear_shape(shade_shape* shape);

// shape properties
void set_points(shade_shape* shape, const vector<int>& points);
void set_lines(shade_shape* shape, const vector<vec2i>& lines);
void set_triangles(shade_shape* shape, const vector<vec3i>& triangles);
void set_quads(shade_shape* shape, const vector<vec4i>& quads);
void set_positions(shade_shape* shape, const vector<vec3f>& positions);
void set_normals(shade_shape* shape, const vector<vec3f>& normals);
void set_texcoords(shade_shape* shape, const vector<vec2f>& texcoords);
void set_colors(shade_shape* shape, const vector<vec4f>& colors);
void set_tangents(shade_shape* shape, const vector<vec4f>& tangents);
void set_instances(
    shade_shape* shape, const vector<vec3f>& froms, const vector<vec3f>& tos);

// get shaoe properties
const ogl_arraybuffer* get_positions(const shade_shape* shape);
const ogl_arraybuffer* get_normals(const shade_shape* shape);
const ogl_arraybuffer* get_texcoords(const shade_shape* shape);
const ogl_arraybuffer* get_colors(const shade_shape* shape);
const ogl_arraybuffer* get_tangents(const shade_shape* shape);

// instance properties
void set_frame(shade_instance* instance, const frame3f& frame);
void set_shape(shade_instance* instance, shade_shape* shape);
void set_material(shade_instance* instance, shade_material* material);
void set_hidden(shade_instance* instance, bool hidden);
void set_highlighted(shade_instance* instance, bool highlighted);

// check if initialized
bool is_initialized(const shade_environment* environment);
// clear environment
void clear_environment(shade_environment* environment);

// environment properties
void set_frame(shade_environment* environment, const frame3f& frame);
void set_emission(shade_environment* environment, const vec3f& emission,
    shade_texture* emission_tex = nullptr);

// shortcuts
shade_camera*   add_camera(shade_scene* scene, const frame3f& frame, float lens,
      float aspect, float film = 0.036, float near = 0.001, float far = 10000);
shade_material* add_material(shade_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    shade_texture* emission_tex = nullptr, shade_texture* color_tex = nullptr,
    shade_texture* specular_tex  = nullptr,
    shade_texture* metallic_tex  = nullptr,
    shade_texture* roughness_tex = nullptr,
    shade_texture* normalmap_tex = nullptr);
shade_shape*    add_shape(shade_scene* scene, const vector<int>& points,
       const vector<vec2i>& lines, const vector<vec3i>& triangles,
       const vector<vec4i>& quads, const vector<vec3f>& positions,
       const vector<vec3f>& normals, const vector<vec2f>& texcoords,
       const vector<vec4f>& colors, bool edges = false);
shade_instance* add_instance(shade_scene* scene, const frame3f& frame,
    shade_shape* shape, shade_material* material, bool hidden = false,
    bool highlighted = false);
shade_environment* add_environment(shade_scene* scene, const frame3f& frame,
    const vec3f& emission, shade_texture* emission_tex = nullptr);

// internal drawing functions
struct shade_view {
  frame3f camera_frame      = {};
  mat4f   view_matrix       = {};
  mat4f   projection_matrix = {};
};

void set_view_uniforms(ogl_program* program, const shade_view& view);
void set_instance_uniforms(ogl_program* program, const frame3f& frame,
    const shade_shape* shape, const shade_material* material,
    const shade_params& params);
void set_camlight_uniforms(
    ogl_program* program, const shade_scene* scene, const shade_view& view);
void set_envlight_uniforms(
    ogl_program* program, const shade_scene* scene, const shade_view& view);

void draw_instances(
    shade_scene* scene, const shade_view& view, const shade_params& params);
void draw_environments(
    shade_scene* scene, const shade_view& view, const shade_params& params);

void draw_scene(shade_scene* scene, shade_camera* camera, const vec4i& viewport,
    const shade_params& params);

// read-only access to defualt shader code
const char* shade_scene_vertex();
const char* shade_camlight_fragment();
const char* shade_envlight_fragment();
const char* shade_enivronment_fragment();

}  // namespace yocto

#endif
