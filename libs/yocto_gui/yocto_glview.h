//
// Yocto/ImageViewer: Simpler image viewer.
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

#ifndef _YOCTO_IMAGEVIEWER_
#define _YOCTO_IMAGEVIEWER_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_trace.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_shade.h>

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void view_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show a set of images
void view_images(const string& title, const vector<string>& names,
    const vector<color_image>& images);

// Open a window and show an image for color grading
void colorgrade_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    const trace_params& params = {}, bool print = true, bool edit = false);

using glview_callback =
    std::function<void(gui_window* win, const gui_input& input,
        vector<int>& updated_shapes, vector<int>& updated_textures)>;

// Open a window and show an scene via OpenGL shading
void glview_scene(const string& title, const string& name, scene_model& scene,
    const shade_params&    params            = {},
    const glview_callback& widgets_callback  = {},
    const glview_callback& uiupdate_callback = {},
    const glview_callback& update_callback   = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL image data
struct glimage_state {
  // image properties
  int width  = 0;
  int height = 0;

  // Opengl state
  uint texture     = 0;  // texture
  uint program     = 0;  // program
  uint vertex      = 0;
  uint fragment    = 0;
  uint vertexarray = 0;  // vertex
  uint positions   = 0;
  uint triangles   = 0;  // elements
};

// create image drawing program
bool init_image(glimage_state& glimage);

// clear image
void clear_image(glimage_state& glimage);

// update image data
void set_image(glimage_state& glimage, const color_image& image);

// OpenGL image drawing params
struct glimage_params {
  vec2i window      = {512, 512};
  vec4i framebuffer = {0, 0, 512, 512};
  vec2f center      = {0, 0};
  float scale       = 1;
  bool  fit         = true;
  bool  checker     = true;
  float border_size = 2;
  vec4f background  = {0.15f, 0.15f, 0.15f, 1.0f};
};

// draw image
void draw_image(glimage_state& image, const glimage_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl caemra
struct glscene_camera {
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
struct glscene_texture : ogl_texture {};

// Opengl material
struct glscene_material {
  // material
  vec3f            emission      = {0, 0, 0};
  vec3f            color         = {0, 0, 0};
  float            metallic      = 0;
  float            roughness     = 0;
  float            specular      = 0;
  float            opacity       = 1;
  gltexture_handle emission_tex  = glinvalid_handle;
  gltexture_handle color_tex     = glinvalid_handle;
  gltexture_handle metallic_tex  = glinvalid_handle;
  gltexture_handle roughness_tex = glinvalid_handle;
  gltexture_handle specular_tex  = glinvalid_handle;
  gltexture_handle opacity_tex   = glinvalid_handle;
  gltexture_handle normal_tex    = glinvalid_handle;
  bool             unlit         = false;
};

// Opengl shape
struct glscene_shape : ogl_shape {};

// Opengl instance
struct glscene_instance {
  // instance properties
  frame3f           frame       = identity3x4f;
  glshape_handle    shape       = glinvalid_handle;
  glmaterial_handle material    = glinvalid_handle;
  bool              hidden      = false;
  bool              highlighted = false;
};

// Opengl environment
struct glscene_environment {
  // environment properties
  frame3f          frame        = identity3x4f;
  vec3f            emission     = {1, 1, 1};
  gltexture_handle emission_tex = glinvalid_handle;

  // drawing data
  glshape_handle   envlight_shape   = glinvalid_handle;
  glcubemap_handle envlight_cubemap = glinvalid_handle;

  // envlight precomputed data
  glcubemap_handle envlight_diffuse_  = glinvalid_handle;
  glcubemap_handle envlight_specular_ = glinvalid_handle;
  gltexture_handle envlight_brdflut_  = glinvalid_handle;
};

// Opengl scene
struct glscene_state {
  // scene objects
  vector<glscene_camera>      cameras      = {};
  vector<glscene_instance>    instances    = {};
  vector<glscene_shape>       shapes       = {};
  vector<glscene_material>    materials    = {};
  vector<glscene_texture>     textures     = {};
  vector<glscene_environment> environments = {};

  // data for envmaps
  vector<glscene_shape> envlight_shapes    = {};
  vector<ogl_cubemap>   envlight_cubemaps  = {};
  vector<ogl_cubemap>   envlight_diffuses  = {};
  vector<ogl_cubemap>   envlight_speculars = {};
  vector<ogl_texture>   envlight_brdfluts  = {};

  // programs
  ogl_program environment_program;
  ogl_program instance_program;

  // disable copy construction
  glscene_state()                     = default;
  glscene_state(const glscene_state&) = delete;
  glscene_state& operator=(const glscene_state&) = delete;

  // cleanup
  ~glscene_state();
};

// Shading type
enum struct glscene_lighting_type { envlight, camlight, eyelight };

// Shading labels
const auto glscene_lighting_labels =
    vector<pair<glscene_lighting_type, string>>{
        {glscene_lighting_type::envlight, "envlight"},
        {glscene_lighting_type::camlight, "camlight"},
        {glscene_lighting_type::eyelight, "eyelight"}};

// Shading name
const auto glscene_lighting_names = vector<string>{
    "envlight", "camlight", "eyelight"};

// Draw options
struct glscene_params {
  int                   camera           = 0;
  int                   resolution       = 1280;
  bool                  wireframe        = false;
  glscene_lighting_type lighting         = glscene_lighting_type::camlight;
  float                 exposure         = 0;
  float                 gamma            = 2.2f;
  bool                  faceted          = false;
  bool                  double_sided     = true;
  bool                  non_rigid_frames = true;
  float                 near             = 0.01f;
  float                 far              = 10000.0f;
  bool                  hide_environment = false;
  vec4f                 background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void init_scene(glscene_state& scene, bool instanced_drawing = false);
bool is_initialized(const glscene_state& scene);

// Initialize data for environment lighting
void init_environments(glscene_state& scene, bool precompute_envlight = true);

// Check if we have an envlight
bool has_envlight(const glscene_state& scene);

// Clear an OpenGL scene
void clear_scene(glscene_state& scene);

// add scene elements
glcamera_handle      add_camera(glscene_state& scene);
gltexture_handle     add_texture(glscene_state& scene);
glmaterial_handle    add_material(glscene_state& scene);
glshape_handle       add_shape(glscene_state& scene);
glinstance_handle    add_instance(glscene_state& scene);
glenvironment_handle add_environment(glscene_state& scene);

// camera properties
void set_frame(glscene_camera& camera, const frame3f& frame);
void set_lens(glscene_camera& camera, float lens, float aspect, float film);
void set_nearfar(glscene_camera& camera, float near, float far);

// material properties
void set_emission(glscene_material& material, const vec3f& emission,
    gltexture_handle emission_tex = glinvalid_handle);
void set_color(glscene_material& material, const vec3f& color,
    gltexture_handle color_tex = glinvalid_handle);
void set_metallic(glscene_material& material, float metallic,
    gltexture_handle metallic_tex = glinvalid_handle);
void set_roughness(glscene_material& material, float roughness,
    gltexture_handle roughness_tex = glinvalid_handle);
void set_specular(glscene_material& material, float specular,
    gltexture_handle specular_tex = glinvalid_handle);
void set_opacity(glscene_material& material, float opacity,
    gltexture_handle opacity_tex = glinvalid_handle);
void set_normalmap(glscene_material& material, gltexture_handle normal_tex);
void set_unlit(glscene_material& material, bool unlit);

// shape properties
void set_points(glscene_shape& shape, const vector<int>& points);
void set_lines(glscene_shape& shape, const vector<vec2i>& lines);
void set_triangles(glscene_shape& shape, const vector<vec3i>& triangles);
void set_quads(glscene_shape& shape, const vector<vec4i>& quads);
void set_positions(glscene_shape& shape, const vector<vec3f>& positions);
void set_normals(glscene_shape& shape, const vector<vec3f>& normals);
void set_texcoords(glscene_shape& shape, const vector<vec2f>& texcoords);
void set_colors(glscene_shape& shape, const vector<vec4f>& colors);
void set_tangents(glscene_shape& shape, const vector<vec4f>& tangents);
void set_instances(
    glscene_shape& shape, const vector<vec3f>& froms, const vector<vec3f>& tos);

// set point size
void set_point_size(glscene_shape& shape, float point_size);

// get shape properties
bool                   has_normals(const glscene_shape& shape);
const ogl_arraybuffer& get_positions(const glscene_shape& shape);
const ogl_arraybuffer& get_normals(const glscene_shape& shape);
const ogl_arraybuffer& get_texcoords(const glscene_shape& shape);
const ogl_arraybuffer& get_colors(const glscene_shape& shape);
const ogl_arraybuffer& get_tangents(const glscene_shape& shape);

// instance properties
void set_frame(glscene_instance& instance, const frame3f& frame);
void set_shape(glscene_instance& instance, glshape_handle shape);
void set_material(glscene_instance& instance, glmaterial_handle material);
void set_hidden(glscene_instance& instance, bool hidden);
void set_highlighted(glscene_instance& instance, bool highlighted);

// environment properties
void set_frame(glscene_environment& environment, const frame3f& frame);
void set_emission(glscene_environment& environment, const vec3f& emission,
    gltexture_handle emission_tex = glinvalid_handle);

// shortcuts
glcamera_handle      add_camera(glscene_state& scene, const frame3f& frame,
         float lens, float aspect, float film = 0.036, float near = 0.001,
         float far = 10000);
glmaterial_handle    add_material(glscene_state& scene, const vec3f& emission,
       const vec3f& color, float specular, float metallic, float roughness,
       gltexture_handle emission_tex  = glinvalid_handle,
       gltexture_handle color_tex     = glinvalid_handle,
       gltexture_handle specular_tex  = glinvalid_handle,
       gltexture_handle metallic_tex  = glinvalid_handle,
       gltexture_handle roughness_tex = glinvalid_handle,
       gltexture_handle normalmap_tex = glinvalid_handle);
glshape_handle       add_shape(glscene_state& scene, const vector<int>& points,
          const vector<vec2i>& lines, const vector<vec3i>& triangles,
          const vector<vec4i>& quads, const vector<vec3f>& positions,
          const vector<vec3f>& normals, const vector<vec2f>& texcoords,
          const vector<vec4f>& colors, bool edges = false);
glinstance_handle    add_instance(glscene_state& scene, const frame3f& frame,
       glmaterial_handle shape, glmaterial_handle material, bool hidden = false,
       bool highlighted = false);
glenvironment_handle add_environment(glscene_state& scene, const frame3f& frame,
    const vec3f& emission, gltexture_handle emission_tex = glinvalid_handle);

// draw scene
void draw_scene(
    glscene_state& scene, const vec4i& viewport, const glscene_params& params);

}  // namespace yocto

#endif
