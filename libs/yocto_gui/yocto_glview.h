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

// Opengl texture
struct glscene_texture {
  // texture properties
  int width  = 0;
  int height = 0;

  // opengl state
  uint texture = 0;
};

// Create texture
void set_texture(glscene_texture& gltexture, const scene_texture& texture);

// Clean texture
void clear_texture(glscene_texture& gltexture);

// Opengl shape
struct glscene_shape {
  // OpenGL objects
  vector<ogl_arraybuffer> vertex_buffers = {};
  ogl_elementbuffer       index_buffer   = {};
  ogl_element_type        elements       = ogl_element_type::triangles;
  size_t                  num_instances  = 0;
  float                   point_size     = 1;

  // OpenGl state
  uint shape_id = 0;

  // Disable copy construction
  glscene_shape()                     = default;
  glscene_shape(const glscene_shape&) = delete;
  glscene_shape& operator=(const glscene_shape&) = delete;
  glscene_shape(glscene_shape&& other) {
    vertex_buffers.swap(other.vertex_buffers);
    std::swap(index_buffer, other.index_buffer);
    std::swap(elements, other.elements);
    std::swap(num_instances, other.num_instances);
    std::swap(point_size, other.point_size);
    std::swap(shape_id, other.shape_id);
  }
};

// Create shape
void set_shape(glscene_shape& glshape, const scene_shape& shape);

// Clean shape
void clear_shape(glscene_shape& glshape);

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
  vector<glscene_shape>       shapes       = {};
  vector<glscene_texture>     textures     = {};
  vector<glscene_environment> environments = {};

  // data for envmaps
  vector<ogl_shape>   envlight_shapes    = {};
  vector<ogl_cubemap> envlight_cubemaps  = {};
  vector<ogl_cubemap> envlight_diffuses  = {};
  vector<ogl_cubemap> envlight_speculars = {};
  vector<ogl_texture> envlight_brdfluts  = {};

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
gltexture_handle     add_texture(glscene_state& scene);
glshape_handle       add_shape(glscene_state& scene);
glenvironment_handle add_environment(glscene_state& scene);

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

// environment properties
void set_frame(glscene_environment& environment, const frame3f& frame);
void set_emission(glscene_environment& environment, const vec3f& emission,
    gltexture_handle emission_tex = glinvalid_handle);

// shortcuts
glshape_handle       add_shape(glscene_state& scene, const vector<int>& points,
          const vector<vec2i>& lines, const vector<vec3i>& triangles,
          const vector<vec4i>& quads, const vector<vec3f>& positions,
          const vector<vec3f>& normals, const vector<vec2f>& texcoords,
          const vector<vec4f>& colors, bool edges = false);
glenvironment_handle add_environment(glscene_state& scene, const frame3f& frame,
    const vec3f& emission, gltexture_handle emission_tex = glinvalid_handle);

// draw scene
void draw_scene(glscene_state& glscene, const scene_model& scene,
    const vec4i& viewport, const glscene_params& params);

}  // namespace yocto

#endif
