//
// Yocto/OpenGL: Utilities for writing native GUIs with OpenGL.
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

#ifndef _YOCTO_OPENGL_
#define _YOCTO_OPENGL_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

#include <array>
#include <string>
#include <vector>

// forward declaration
struct GLFWwindow;

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
// LOW-LEVEL OPENGL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Commands to setup the opengl context and issue gpu operations.
bool init_ogl(string& error);
void assert_ogl_error();
bool check_ogl_error(string& error);
void clear_ogl_framebuffer(const vec4f& color, bool clear_depth = true);
void set_ogl_viewport(const vec4i& viewport);
void set_ogl_viewport(const vec2i& viewport);
void set_ogl_wireframe(bool enabled);
void set_ogl_blending(bool enabled);
void set_ogl_point_size(int size);

// OpenGL texture
struct ogl_texture {
  // Texture properties
  vec2i size      = {0, 0};
  int   nchannels = 0;
  bool  is_srgb   = false;
  bool  is_float  = false;
  bool  linear    = false;
  bool  mipmap    = false;

  // OpenGL state
  uint texture_id = 0;

  // ensuring no copies
  ogl_texture() {}
  ogl_texture(const ogl_texture&) = delete;
  ogl_texture& operator=(ogl_texture&) = delete;
};

// set texture
void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const byte* img, bool as_srgb = false, bool linear = true,
    bool mipmap = true);
void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const float* img, bool as_float = false, bool linear = true,
    bool mipmap = true);

// check if texture is initialized
bool is_initialized(ogl_texture* texture);

// clear texture
void clear_texture(ogl_texture* texture);

// set texture
void set_texture(ogl_texture* texture, const image<vec4b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec4f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec3b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec3f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<byte>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<float>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);

struct ogl_cubemap {
  // Cubemap properties
  int  size      = 0;
  int  nchannels = 0;
  bool is_srgb   = false;
  bool is_float  = false;
  bool linear    = false;
  bool mipmap    = false;

  // OpenGL state
  uint cubemap_id = 0;
};

// set cubemap
void set_cubemap(ogl_cubemap* cubemap, int size, int nchannels,
    const array<byte*, 6>& img, bool as_srgb = false, bool linear = true,
    bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, int size, int nchannels,
    const array<float*, 6>& img, bool as_float = false, bool linear = true,
    bool mipmap = true);

template <typename T>
void set_cubemap(ogl_cubemap* cubemap, int size, int nchannels,
    bool as_float = false, bool linear = true, bool mipmap = true) {
  auto img = array<T*, 6>{0, 0, 0, 0, 0, 0};
  set_cubemap(cubemap, size, nchannels, img, as_float, mipmap);
}

// check if cubemap is initialized
bool is_initialized(ogl_cubemap* cubemap);

// clear cubemap
void clear_cubemap(ogl_cubemap* cubemap);

void set_cubemap(ogl_cubemap* cubemap, const array<image<vec4b>, 6>& img,
    int nchannels, bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec4f>, 6>& img,
    int nchannels, bool as_float = false, bool linear = true,
    bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec3b>, 6>& img,
    int nchannels, bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec3f>, 6>& img,
    int nchannels, bool as_float = false, bool linear = true,
    bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, const array<image<byte>, 6>& img,
    int nchannels, bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_cubemap(ogl_cubemap* cubemap, const array<image<float>, 6>& img,
    int nchannels, bool as_float = false, bool linear = true,
    bool mipmap = true);

// Opengl array/element buffer
struct ogl_arraybuffer {
  // buffer data
  size_t size    = 0;
  int    esize   = 0;
  bool   dynamic = false;
  // OpenGL state
  uint buffer_id = 0;
};

// set buffer
void set_arraybuffer(ogl_arraybuffer* buffer, size_t size, int esize,
    const float* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(ogl_arraybuffer* buffer);

// clear buffer
void clear_arraybuffer(ogl_arraybuffer* buffer);

// set buffer
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<float>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec2f>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec3f>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec4f>& data, bool dynamic = false);

// Opengl draw elements
enum struct ogl_element_type { points, lines, triangles };

// Opengl array/element buffer
struct ogl_elementbuffer {
  // buffer data
  size_t           size    = 0;
  ogl_element_type element = ogl_element_type::points;
  bool             dynamic = false;
  // OpenGL state
  uint buffer_id = 0;
};

// set buffer
void set_elementbuffer(ogl_elementbuffer* buffer, size_t size,
    ogl_element_type element, const int* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(ogl_elementbuffer* buffer);

// clear buffer
void clear_elementbuffer(ogl_elementbuffer* buffer);

// set buffer
void set_elementbuffer(
    ogl_elementbuffer* buffer, const vector<int>& points, bool dynamic = false);
void set_elementbuffer(ogl_elementbuffer* buffer, const vector<vec2i>& lines,
    bool dynamic = false);
void set_elementbuffer(ogl_elementbuffer* buffer,
    const vector<vec3i>& triangles, bool dynamic = false);

// Opengl program
struct ogl_program {
  // program code
  string vertex_code;
  string fragment_code;
  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
};

// initialize program
bool init_program(ogl_program* program, const string& vertex,
    const string& fragment, string& error, string& errorlog);
bool is_initialized(const ogl_program* program);

// clear program
void clear_program(ogl_program* program);

// bind program
void bind_program(ogl_program* program);
// unbind program
void unbind_program(ogl_program* program);
// unbind program
void unbind_program();

// get uniform location
int get_uniform_location(ogl_program* program, const char* name);

// set uniforms
void set_uniform(ogl_program* program, int location, int value);
void set_uniform(ogl_program* program, int location, const vec2i& value);
void set_uniform(ogl_program* program, int location, const vec3i& value);
void set_uniform(ogl_program* program, int location, const vec4i& value);
void set_uniform(ogl_program* program, int location, float value);
void set_uniform(ogl_program* program, int location, const vec2f& value);
void set_uniform(ogl_program* program, int location, const vec3f& value);
void set_uniform(ogl_program* program, int location, const vec4f& value);
void set_uniform(ogl_program* program, int location, const mat2f& value);
void set_uniform(ogl_program* program, int location, const mat3f& value);
void set_uniform(ogl_program* program, int location, const mat4f& value);
void set_uniform(ogl_program* program, int location, const frame2f& value);
void set_uniform(ogl_program* program, int location, const frame3f& value);
template <typename T>
inline void set_uniform(
    ogl_program* program, const char* name, const T& value) {
  return set_uniform(program, get_uniform_location(program, name), value);
}

// set uniform texture
void set_uniform(
    ogl_program* program, int location, const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, const char* name,
    const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, int location, int location_on,
    const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, const char* name, const char* name_on,
    const ogl_texture* texture, int unit);

// set uniform cubemap
void set_uniform(
    ogl_program* program, int location, const ogl_cubemap* cubemap, int unit);
void set_uniform(ogl_program* program, const char* name,
    const ogl_cubemap* cubemap, int unit);
void set_uniform(ogl_program* program, int location, int location_on,
    const ogl_cubemap* cubemap, int unit);
void set_uniform(ogl_program* program, const char* name, const char* name_on,
    const ogl_cubemap* cubemap, int unit);

// get attribute location
int get_attribute_location(ogl_program* program, const char* name);

// set vertex attributes
void set_attribute(ogl_program* program, int location, ogl_arraybuffer* buffer);
void set_attribute(
    ogl_program* program, const char* name, ogl_arraybuffer* buffer);

// set vertex attributes
void set_attribute(ogl_program* program, int location, float value);
void set_attribute(ogl_program* program, int location, const vec2f& value);
void set_attribute(ogl_program* program, int location, const vec3f& value);
void set_attribute(ogl_program* program, int location, const vec4f& value);
template <typename T>
inline void set_attribute(
    ogl_program* program, const char* name, const T& value) {
  return set_attribute(program, get_attribute_location(program, name), value);
}

// set vertex attributes
template <typename T>
inline void set_attribute(ogl_program* program, int location,
    ogl_arraybuffer* buffer, const T& value) {
  if (buffer && is_initialized(buffer)) {
    return set_attribute(program, location, buffer);
  } else {
    set_attribute(program, location, value);
  }
}
template <typename T>
inline void set_attribute(ogl_program* program, const char* name,
    ogl_arraybuffer* buffer, const T& def) {
  set_attribute(program, get_attribute_location(program, name), buffer, def);
}

// draw elements
void draw_elements(ogl_elementbuffer* buffer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FRAMEBUFFERS
// -----------------------------------------------------------------------------
namespace yocto {

struct ogl_framebuffer {
  vec2i size            = {0, 0};
  uint  framebuffer_id  = 0;
  uint  renderbuffer_id = 0;

  static inline uint bound_framebuffer_id = 0;
  // TODO(giacomo): Maybe the following is better
  // static inline ogl_framebuffer* bound_framebuffer = nullptr;
};

void set_framebuffer(ogl_framebuffer* framebuffer, const vec2i& size);

void set_framebuffer_texture(const ogl_framebuffer* framebuffer,
    const ogl_texture* texture, uint mipmap_level = 0);

void set_framebuffer_texture(const ogl_framebuffer* framebuffer,
    const ogl_cubemap* cubemap, uint face, uint mipmap_level = 0);

void bind_framebuffer(const ogl_framebuffer* target);
void unbind_framebuffer();
void clear_framebuffer(ogl_framebuffer* target);
}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPES
// -----------------------------------------------------------------------------
namespace yocto {
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

  uint shape_id = 0;

  ogl_shape() {}
  ogl_shape(const ogl_shape&) = delete;
  ogl_shape& operator=(const ogl_shape&) = delete;
  ~ogl_shape();
};

void set_shape(ogl_shape* shape);
void clear_shape(ogl_shape* shape);
void bind_shape(const ogl_shape* shape);

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
void set_colors(ogl_shape* shape, const vector<vec4f>& colors);
void set_tangents(ogl_shape* shape, const vector<vec4f>& tangents);

void draw_shape(const ogl_shape* shape);

ogl_shape* cube_shape();
ogl_shape* quad_shape();

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL image data
struct ogl_image {
  ogl_image() {}
  ogl_image(const ogl_image&) = delete;
  ogl_image& operator=(const ogl_image&) = delete;
  ~ogl_image();

  ogl_program* program = new ogl_program{};
  ogl_texture* texture = new ogl_texture{};
  ogl_shape*   quad    = new ogl_shape{};
};

// create image drawing program
bool init_image(ogl_image* oimg);
bool is_initialized(const ogl_image* oimg);

// clear image
void clear_image(ogl_image* oimg);

// update image data
void set_image(ogl_image* oimg, const image<vec4f>& img, bool linear = false,
    bool mipmap = false);
void set_image(ogl_image* oimg, const image<vec4b>& img, bool linear = false,
    bool mipmap = false);

// OpenGL image drawing params
struct ogl_image_params {
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
void draw_image(ogl_image* image, const ogl_image_params& params);

}  // namespace yocto

#endif
