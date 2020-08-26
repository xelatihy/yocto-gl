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

struct gui_texture : ogl_texture {};

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

struct gui_shape : ogl_shape {};

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
ogl_arraybuffer* get_positions(gui_shape* shape);
ogl_arraybuffer* get_normals(gui_shape* shape);
ogl_arraybuffer* get_texcoords(gui_shape* shape);
ogl_arraybuffer* get_colors(gui_shape* shape);
ogl_arraybuffer* get_tangents(gui_shape* shape);

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

// Opengl scene
struct gui_scene {
  // scene objects
  vector<gui_camera*>   cameras   = {};
  vector<gui_instance*> instances = {};
  vector<gui_shape*>    shapes    = {};
  vector<gui_material*> materials = {};
  vector<gui_texture*>  textures  = {};

  // programs
  ogl_program* eyelight_program    = new ogl_program{};
  ogl_program* ibl_program         = new ogl_program{};
  ogl_program* environment_program = new ogl_program{};

  // IBL baked data
  ogl_cubemap* environment_cubemap = new ogl_cubemap{};
  ogl_cubemap* diffuse_cubemap     = new ogl_cubemap{};
  ogl_cubemap* specular_cubemap    = new ogl_cubemap{};
  gui_texture* brdf_lut            = new gui_texture{};

  ~gui_scene();
};

// Shading type
enum struct gui_lighting_type {
  environment,
  eyelight,
  // scene_lights
};

// Shading name
const auto gui_lighting_names = vector<string>{"environment", "camera_lights"};

// Draw options
struct gui_scene_params {
  int               resolution       = 1280;
  bool              wireframe        = false;
  gui_lighting_type lighting         = gui_lighting_type::eyelight;
  float             exposure         = 0;
  float             gamma            = 2.2f;
  bool              faceted          = false;
  bool              double_sided     = true;
  bool              non_rigid_frames = true;
  float             near             = 0.01f;
  float             far              = 10000.0f;
  vec4f             background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

struct gui_scene_view {
  frame3f          camera_frame      = {};
  mat4f            view_matrix       = {};
  mat4f            projection_matrix = {};
  gui_scene_params params            = {};
};

// Initialize an OpenGL scene
void init_scene(gui_scene* scene);
bool is_initialized(const gui_scene* scene);

// Initialize data for image based lighting
void init_ibl_data(
    gui_scene* scene, const gui_texture* environment, const vec3f& emission);

// Clear an OpenGL scene
void clear_scene(gui_scene* scene);

// add scene elements
gui_camera*   add_camera(gui_scene* scene);
gui_texture*  add_texture(gui_scene* scene);
gui_material* add_material(gui_scene* scene);
gui_shape*    add_shape(gui_scene* scene);
gui_instance* add_instance(gui_scene* scene);

// camera properties
void set_frame(gui_camera* camera, const frame3f& frame);
void set_lens(gui_camera* camera, float lens, float aspect, float film);
void set_nearfar(gui_camera* camera, float near, float far);

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

gui_shape* add_shape(gui_scene* scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec4f>& colors, bool edges = false);

// instance properties
void set_frame(gui_instance* instance, const frame3f& frame);
void set_shape(gui_instance* instance, gui_shape* shape);
void set_material(gui_instance* instance, gui_material* material);
void set_hidden(gui_instance* instance, bool hidden);
void set_highlighted(gui_instance* instance, bool highlighted);

// shortcuts
gui_camera*   add_camera(gui_scene* scene, const frame3f& frame, float lens,
      float aspect, float film = 0.036, float near = 0.001, float far = 10000);
gui_material* add_material(gui_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    gui_texture* emission_tex = nullptr, gui_texture* color_tex = nullptr,
    gui_texture* specular_tex = nullptr, gui_texture* metallic_tex = nullptr,
    gui_texture* roughness_tex = nullptr, gui_texture* normalmap_tex = nullptr);
gui_instance* add_instance(gui_scene* scene, const frame3f& frame,
    gui_shape* shape, gui_material* material, bool hidden = false,
    bool highlighted = false);

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
