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

// forward declaration
struct GLFWwindow;

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
struct ogl_camera {
  frame3f frame    = identity3x4f;
  float   lens     = 0.050;
  float   aspect   = 1.000;
  float   film     = 0.036;
  float   near     = 0.001;
  float   far      = 10000;
  float   aperture = 0;
  float   focus    = 0;
};

// Opengl material
struct ogl_material {
  // material
  vec3f        emission      = {0, 0, 0};
  vec3f        color         = {0, 0, 0};
  float        metallic      = 0;
  float        roughness     = 0;
  float        specular      = 0;
  float        opacity       = 1;
  ogl_texture* emission_tex  = nullptr;
  ogl_texture* color_tex     = nullptr;
  ogl_texture* metallic_tex  = nullptr;
  ogl_texture* roughness_tex = nullptr;
  ogl_texture* specular_tex  = nullptr;
  ogl_texture* opacity_tex   = nullptr;
  ogl_texture* normal_tex    = nullptr;
};

// Opengl shape
struct ogl_shape {
  // vertex buffers
  ogl_arraybuffer*   positions      = new ogl_arraybuffer{};
  ogl_arraybuffer*   normals        = new ogl_arraybuffer{};
  ogl_arraybuffer*   texcoords      = new ogl_arraybuffer{};
  ogl_arraybuffer*   colors         = new ogl_arraybuffer{};
  ogl_arraybuffer*   tangents       = new ogl_arraybuffer{};
  ogl_elementbuffer* points         = new ogl_elementbuffer{};
  ogl_elementbuffer* lines          = new ogl_elementbuffer{};
  ogl_elementbuffer* triangles      = new ogl_elementbuffer{};
  ogl_elementbuffer* quads          = new ogl_elementbuffer{};
  ogl_elementbuffer* edges          = new ogl_elementbuffer{};
  float              points_size    = 10;
  float              line_thickness = 4;

  ogl_shape() {}
  ogl_shape(const ogl_shape&) = delete;
  ogl_shape& operator=(const ogl_shape&) = delete;
  ~ogl_shape();
};

// Opengl instance
struct ogl_instance {
  // instance properties
  frame3f       frame       = identity3x4f;
  ogl_shape*    shape       = nullptr;
  ogl_material* material    = nullptr;
  bool          hidden      = false;
  bool          highlighted = false;
};

// Light type
enum struct ogl_light_type { point = 0, directional };

// Opengl light
struct ogl_light {
  vec3f          position = {0, 0, 0};
  vec3f          emission = {0, 0, 0};
  ogl_light_type type     = ogl_light_type::point;
  bool           camera   = false;
};

// Opengl scene
struct ogl_scene {
  ogl_scene() {}
  ogl_scene(const ogl_scene&) = delete;
  ogl_scene& operator=(const ogl_scene&) = delete;
  ~ogl_scene();

  // scene objects
  vector<ogl_camera*>   cameras   = {};
  vector<ogl_instance*> instances = {};
  vector<ogl_shape*>    shapes    = {};
  vector<ogl_material*> materials = {};
  vector<ogl_texture*>  textures  = {};
  vector<ogl_light*>    lights    = {};

  // OpenGL state
  ogl_program* program = new ogl_program{};

  ogl_program* environment_program = new ogl_program{};
  ogl_cubemap* environment_cubemap = new ogl_cubemap{};
  ogl_cubemap* irradiance_map      = new ogl_cubemap{};
  ogl_cubemap* prefiltered_map     = new ogl_cubemap{};
  ogl_texture* brdf_lut            = new ogl_texture{};
};

// Shading type
enum struct ogl_shading_type { lights, eyelight, camlights };

// Shading name
const auto ogl_shading_names = vector<string>{
    "lights", "eyelight", "camlights"};

// Draw options
struct ogl_scene_params {
  int              resolution       = 1280;
  bool             wireframe        = false;
  bool             edges            = false;
  float            edge_offset      = 0.01f;
  ogl_shading_type shading          = ogl_shading_type::camlights;
  float            exposure         = 0;
  float            gamma            = 2.2f;
  vec3f            ambient          = {0, 0, 0};
  bool             double_sided     = true;
  bool             non_rigid_frames = true;
  float            near             = 0.01f;
  float            far              = 10000.0f;
  vec4f            background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void init_scene(ogl_scene* scene);
bool is_initialized(const ogl_scene* scene);

// Clear an OpenGL scene
void clear_scene(ogl_scene* scene);

// add scene elements
ogl_camera*   add_camera(ogl_scene* scene);
ogl_texture*  add_texture(ogl_scene* scene);
ogl_material* add_material(ogl_scene* scene);
ogl_shape*    add_shape(ogl_scene* scene);
ogl_instance* add_instance(ogl_scene* scene);
ogl_light*    add_light(ogl_scene* scene);

// camera properties
void set_frame(ogl_camera* camera, const frame3f& frame);
void set_lens(ogl_camera* camera, float lens, float aspect, float film);
void set_nearfar(ogl_camera* camera, float near, float far);

// material properties
void set_emission(ogl_material* material, const vec3f& emission,
    ogl_texture* emission_tex = nullptr);
void set_color(ogl_material* material, const vec3f& color,
    ogl_texture* color_tex = nullptr);
void set_metallic(ogl_material* material, float metallic,
    ogl_texture* metallic_tex = nullptr);
void set_roughness(ogl_material* material, float roughness,
    ogl_texture* roughness_tex = nullptr);
void set_specular(ogl_material* material, float specular,
    ogl_texture* specular_tex = nullptr);
void set_opacity(
    ogl_material* material, float opacity, ogl_texture* opacity_tex = nullptr);
void set_normalmap(ogl_material* material, ogl_texture* normal_tex);

// shape properties
void set_points(ogl_shape* shape, const vector<int>& points);
void set_lines(ogl_shape* shape, const vector<vec2i>& lines);
void set_triangles(ogl_shape* shape, const vector<vec3i>& triangles);
void set_quads(ogl_shape* shape, const vector<vec4i>& quads);
void set_edges(ogl_shape* shape, const vector<vec3i>& triangles,
    const vector<vec4i>& quads);
void set_positions(ogl_shape* shape, const vector<vec3f>& positions);
void set_normals(ogl_shape* shape, const vector<vec3f>& normals);
void set_texcoords(ogl_shape* shape, const vector<vec2f>& texcoords);
void set_colors(ogl_shape* shape, const vector<vec3f>& colors);
void set_tangents(ogl_shape* shape, const vector<vec4f>& tangents);

// instance properties
void set_frame(ogl_instance* instance, const frame3f& frame);
void set_shape(ogl_instance* instance, ogl_shape* shape);
void set_material(ogl_instance* instance, ogl_material* material);
void set_hidden(ogl_instance* instance, bool hidden);
void set_highlighted(ogl_instance* instance, bool highlighted);

// shortcuts
ogl_camera*   add_camera(ogl_scene* scene, const frame3f& frame, float lens,
      float aspect, float film = 0.036, float near = 0.001, float far = 10000);
ogl_material* add_material(ogl_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    ogl_texture* emission_tex = nullptr, ogl_texture* color_tex = nullptr,
    ogl_texture* specular_tex = nullptr, ogl_texture* metallic_tex = nullptr,
    ogl_texture* roughness_tex = nullptr, ogl_texture* normalmap_tex = nullptr);
ogl_shape*    add_shape(ogl_scene* scene, const vector<int>& points,
       const vector<vec2i>& lines, const vector<vec3i>& triangles,
       const vector<vec4i>& quads, const vector<vec3f>& positions,
       const vector<vec3f>& normals, const vector<vec2f>& texcoords,
       const vector<vec3f>& colors, bool edges = false);
ogl_instance* add_instance(ogl_scene* scene, const frame3f& frame,
    ogl_shape* shape, ogl_material* material, bool hidden = false,
    bool highlighted = false);

// light properties
void add_default_lights(ogl_scene* scene);
void set_light(ogl_light* light, const vec3f& position, const vec3f& emission,
    ogl_light_type type, bool camera);

// light size
void clear_lights(ogl_scene* scene);
bool has_max_lights(ogl_scene* scene);

// Draw an OpenGL scene
void draw_scene(ogl_scene* scene, ogl_camera* camera, const vec4i& viewport,
    const ogl_scene_params& params);

}  // namespace yocto

#endif
