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

// Opengl material
struct gui_material {
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

// Opengl instance
struct gui_instance {
  // instance properties
  frame3f       frame       = identity3x4f;
  ogl_shape*    shape       = nullptr;
  gui_material* material    = nullptr;
  bool          hidden      = false;
  bool          highlighted = false;
};

// Light type
enum struct ogl_light_type { point = 0, directional };

// Opengl light
struct gui_light {
  vec3f          position = {0, 0, 0};
  vec3f          emission = {0, 0, 0};
  ogl_light_type type     = ogl_light_type::point;
  bool           camera   = false;
};

// Opengl scene
struct gui_scene {
  gui_scene() {}
  gui_scene(const gui_scene&) = delete;
  gui_scene& operator=(const gui_scene&) = delete;
  ~gui_scene();

  // scene objects
  vector<gui_camera*>   cameras   = {};
  vector<gui_instance*> instances = {};
  vector<ogl_shape*>    shapes    = {};
  vector<gui_material*> materials = {};
  vector<ogl_texture*>  textures  = {};
  vector<gui_light*>    lights    = {};

  // OpenGL state
  ogl_program* lights_program      = new ogl_program{};
  ogl_program* ibl_program         = new ogl_program{};
  ogl_program* environment_program = new ogl_program{};

  // IBL data
  ogl_cubemap* environment_cubemap = new ogl_cubemap{};
  ogl_cubemap* diffuse_cubemap     = new ogl_cubemap{};
  ogl_cubemap* specular_cubemap    = new ogl_cubemap{};
  ogl_texture* brdf_lut            = new ogl_texture{};
};

// Shading type
enum struct gui_shading_type {
  environment,
  eyelight,
  // scene_lights
};

// Shading name
const auto gui_shading_names = vector<string>{"environment", "camera_lights"};

// Draw options
struct gui_scene_params {
  int              resolution       = 1280;
  bool             wireframe        = false;
  bool             edges            = false;
  float            edge_offset      = 0.01f;
  gui_shading_type shading          = gui_shading_type::eyelight;
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
void init_scene(gui_scene* scene);
bool is_initialized(const gui_scene* scene);

namespace ibl {
// Initialize data for image based lighting
void init_ibl_data(gui_scene* scene, const ogl_texture* environment);
}  // namespace ibl

// Clear an OpenGL scene
void clear_scene(gui_scene* scene);

// add scene elements
gui_camera*   add_camera(gui_scene* scene);
ogl_texture*  add_texture(gui_scene* scene);
gui_material* add_material(gui_scene* scene);
ogl_shape*    add_shape(gui_scene* scene);
gui_instance* add_instance(gui_scene* scene);
gui_light*    add_light(gui_scene* scene);

// camera properties
void set_frame(gui_camera* camera, const frame3f& frame);
void set_lens(gui_camera* camera, float lens, float aspect, float film);
void set_nearfar(gui_camera* camera, float near, float far);

// material properties
void set_emission(gui_material* material, const vec3f& emission,
    ogl_texture* emission_tex = nullptr);
void set_color(gui_material* material, const vec3f& color,
    ogl_texture* color_tex = nullptr);
void set_metallic(gui_material* material, float metallic,
    ogl_texture* metallic_tex = nullptr);
void set_roughness(gui_material* material, float roughness,
    ogl_texture* roughness_tex = nullptr);
void set_specular(gui_material* material, float specular,
    ogl_texture* specular_tex = nullptr);
void set_opacity(
    gui_material* material, float opacity, ogl_texture* opacity_tex = nullptr);
void set_normalmap(gui_material* material, ogl_texture* normal_tex);

ogl_shape* add_shape(gui_scene* scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec4f>& colors, bool edges = false);

// instance properties
void set_frame(gui_instance* instance, const frame3f& frame);
void set_shape(gui_instance* instance, ogl_shape* shape);
void set_material(gui_instance* instance, gui_material* material);
void set_hidden(gui_instance* instance, bool hidden);
void set_highlighted(gui_instance* instance, bool highlighted);

// shortcuts
gui_camera*   add_camera(gui_scene* scene, const frame3f& frame, float lens,
      float aspect, float film = 0.036, float near = 0.001, float far = 10000);
gui_material* add_material(gui_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    ogl_texture* emission_tex = nullptr, ogl_texture* color_tex = nullptr,
    ogl_texture* specular_tex = nullptr, ogl_texture* metallic_tex = nullptr,
    ogl_texture* roughness_tex = nullptr, ogl_texture* normalmap_tex = nullptr);
// ogl_shape*    _add_shape(gui_scene* scene, const vector<int>& points,
//        const vector<vec2i>& lines, const vector<vec3i>& triangles,
//        const vector<vec4i>& quads, const vector<vec3f>& positions,
//        const vector<vec3f>& normals, const vector<vec2f>& texcoords,
//        const vector<vec3f>& colors, bool edges = false);
gui_instance* add_instance(gui_scene* scene, const frame3f& frame,
    ogl_shape* shape, gui_material* material, bool hidden = false,
    bool highlighted = false);

// light properties
void add_default_lights(gui_scene* scene);
void set_light(gui_light* light, const vec3f& position, const vec3f& emission,
    ogl_light_type type, bool camera);

// light size
void clear_lights(gui_scene* scene);
bool has_max_lights(gui_scene* scene);

// Draw an OpenGL scene
void draw_scene(gui_scene* scene, gui_camera* camera, const vec4i& viewport,
    const gui_scene_params& params);

}  // namespace yocto

#endif
