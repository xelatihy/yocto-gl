//
// Utilities to use OpenGL 3, GLFW and ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_opengl.h"

#include <algorithm>
#include <cstdarg>
#include <deque>
#include <mutex>

#include "ext/glad/glad.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"

#define CUTE_FILES_IMPLEMENTATION
#include "ext/cute_files.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// OPENGL UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

static void init_glprogram(uint& program_id, uint& vertex_id, uint& fragment_id,
    uint& array_id, const char* vertex, const char* fragment) {
  assert(glGetError() == GL_NO_ERROR);
  glGenVertexArrays(1, &array_id);
  glBindVertexArray(array_id);
  assert(glGetError() == GL_NO_ERROR);

  int  errflags;
  char errbuf[10000];

  // create vertex
  assert(glGetError() == GL_NO_ERROR);
  vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_id, 1, &vertex, NULL);
  glCompileShader(vertex_id);
  glGetShaderiv(vertex_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(vertex_id, 10000, 0, errbuf);
    throw std::runtime_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create fragment
  assert(glGetError() == GL_NO_ERROR);
  fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_id, 1, &fragment, NULL);
  glCompileShader(fragment_id);
  glGetShaderiv(fragment_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(fragment_id, 10000, 0, errbuf);
    throw std::runtime_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create program
  assert(glGetError() == GL_NO_ERROR);
  program_id = glCreateProgram();
  glAttachShader(program_id, vertex_id);
  glAttachShader(program_id, fragment_id);
  glLinkProgram(program_id);
  glValidateProgram(program_id);
  glGetProgramiv(program_id, GL_LINK_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program_id, 10000, 0, errbuf);
    throw std::runtime_error("program not linked with error\n"s + errbuf);
  }
  glGetProgramiv(program_id, GL_VALIDATE_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program_id, 10000, 0, errbuf);
    throw std::runtime_error("program not linked with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void init_glbuffer(
    uint& buffer_id, bool element, int size, int count, const float* array) {
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer_id);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER,
      count * size * sizeof(float), array, GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void init_glbuffer(
    uint& buffer_id, bool element, int size, int count, const int* array) {
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer_id);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER,
      count * size * sizeof(int), array, GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void update_glbuffer(
    uint& buffer_id, bool element, int size, int count, const int* array) {
  assert(glGetError() == GL_NO_ERROR);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferSubData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, 0,
      size * count * sizeof(int), array);
  assert(glGetError() == GL_NO_ERROR);
}

void update_glbuffer(
    uint& buffer_id, bool element, int size, int count, const float* array) {
  assert(glGetError() == GL_NO_ERROR);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferSubData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, 0,
      size * count * sizeof(float), array);
  assert(glGetError() == GL_NO_ERROR);
}

static void update_gltexture(uint& texture_id, const vec2i& size, int nchan,
    const float* img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
      nchan == 4 ? GL_RGBA : GL_RGB, GL_FLOAT, img);
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

static void update_gltexture(uint& texture_id, const vec2i& size, int nchan,
    const byte* img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
      nchan == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, img);
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

static void init_gltexture(uint& texture_id, const vec2i& size, int nchan,
    const float* img, bool as_float, bool linear, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glGenTextures(1, &texture_id);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  if (as_float) {
    glTexImage2D(GL_TEXTURE_2D, 0, nchan == 4 ? GL_RGBA32F : GL_RGB32F, size.x,
        size.y, 0, nchan == 4 ? GL_RGBA : GL_RGB, GL_FLOAT, img);
  } else {
    glTexImage2D(GL_TEXTURE_2D, 0, nchan == 4 ? GL_RGBA : GL_RGB, size.x,
        size.y, 0, nchan == 4 ? GL_RGBA : GL_RGB, GL_FLOAT, img);
  }
  if (mipmap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  }
  assert(glGetError() == GL_NO_ERROR);
}

static void init_gltexture(uint& texture_id, const vec2i& size, int nchan,
    const byte* img, bool as_srgb, bool linear, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glGenTextures(1, &texture_id);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  if (as_srgb) {
    glTexImage2D(GL_TEXTURE_2D, 0, nchan == 4 ? GL_SRGB_ALPHA : GL_SRGB, size.x,
        size.y, 0, nchan == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, img);
  } else {
    glTexImage2D(GL_TEXTURE_2D, 0, nchan == 4 ? GL_RGBA : GL_RGB, size.x,
        size.y, 0, nchan == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, img);
  }
  if (mipmap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  }
  assert(glGetError() == GL_NO_ERROR);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

auto glimage_vertex =
    R"(
#version 330
in vec2 texcoord;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (texcoord - vec2(0.5,0.5)) * image_size * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
    frag_texcoord = texcoord;
}
)";
#if 0
  auto glimage_vertex = R"(
#version 330
in vec2 texcoord;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size, border_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (texcoord - vec2(0.5,0.5)) * (image_size + border_size*2) * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
    frag_texcoord = texcoord;
}
)";
#endif
auto glimage_fragment =
    R"(
#version 330
in vec2 frag_texcoord;
out vec4 frag_color;
uniform sampler2D txt;
void main() {
    frag_color = texture(txt, frag_texcoord);
}
)";
#if 0
auto glimage_fragment = R"(
#version 330
in vec2 frag_texcoord;
out vec4 frag_color;
uniform vec2 image_size, border_size;
uniform float image_scale;
void main() {
    ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
    ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
    ivec2 tile = tilecoord / 16;
    if(imcoord.x <= 0 || imcoord.y <= 0 || 
        imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
    else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
    else frag_color = vec4(0.3,0.3,0.3,1);
}
)";
#endif

unique_ptr<opengl_image> make_glimage() {
  auto glimage = make_unique<opengl_image>();
  init_glimage(glimage.get());
  return glimage;
}
bool is_initialized(const opengl_image* glimage) {
  return (bool)glimage->program_id;
}

opengl_image::~opengl_image() {
  if (program_id) glDeleteProgram(program_id);
  if (vertex_id) glDeleteShader(vertex_id);
  if (fragment_id) glDeleteShader(fragment_id);
  if (array_id) glDeleteVertexArrays(1, &array_id);
  if (texcoords_id) glDeleteBuffers(1, &texcoords_id);
  if (triangles_id) glDeleteBuffers(1, &triangles_id);
  if (texture_id) glDeleteTextures(1, &texture_id);
}

// init image program
void init_glimage(opengl_image* glimage) {
  if (glimage->program_id) return;

  auto texcoords = vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  auto triangles = vector<vec3i>{{0, 1, 2}, {0, 2, 3}};

  init_glprogram(glimage->program_id, glimage->vertex_id, glimage->fragment_id,
      glimage->array_id, glimage_vertex, glimage_fragment);
  glGenBuffers(1, &glimage->texcoords_id);
  glBindBuffer(GL_ARRAY_BUFFER, glimage->texcoords_id);
  glBufferData(GL_ARRAY_BUFFER, texcoords.size() * 2 * sizeof(float),
      texcoords.data(), GL_STATIC_DRAW);
  glGenBuffers(1, &glimage->triangles_id);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glimage->triangles_id);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * 3 * sizeof(int),
      triangles.data(), GL_STATIC_DRAW);
}

// update image data
void set_glimage(
    opengl_image* glimage, const image<vec4f>& img, bool linear, bool mipmap) {
  if (!glimage->texture_id) {
    init_gltexture(glimage->texture_id, img.size(), 4, &img.data()->x, false,
        linear, mipmap);
  } else if (glimage->texture_size != img.size() ||
             glimage->texture_linear != linear ||
             glimage->texture_mipmap != mipmap) {
    glDeleteTextures(1, &glimage->texture_id);
    init_gltexture(glimage->texture_id, img.size(), 4, &img.data()->x, false,
        linear, mipmap);
  } else {
    update_gltexture(
        glimage->texture_id, img.size(), 4, &img.data()->x, mipmap);
  }
  glimage->texture_size   = img.size();
  glimage->texture_linear = linear;
  glimage->texture_mipmap = mipmap;
}
void set_glimage(
    opengl_image* glimage, const image<vec4b>& img, bool linear, bool mipmap) {
  if (!glimage->texture_id) {
    init_gltexture(glimage->texture_id, img.size(), 4, &img.data()->x, false,
        linear, mipmap);
  } else if (glimage->texture_size != img.size() ||
             glimage->texture_linear != linear ||
             glimage->texture_mipmap != mipmap) {
    glDeleteTextures(1, &glimage->texture_id);
    init_gltexture(glimage->texture_id, img.size(), 4, &img.data()->x, false,
        linear, mipmap);
  } else {
    update_gltexture(
        glimage->texture_id, img.size(), 4, &img.data()->x, mipmap);
  }
  glimage->texture_size   = img.size();
  glimage->texture_linear = linear;
  glimage->texture_mipmap = mipmap;
}

// draw image
void draw_glimage(opengl_image* glimage, const draw_glimage_params& params) {
  assert(glGetError() == GL_NO_ERROR);
  glViewport(params.framebuffer.x, params.framebuffer.y, params.framebuffer.z,
      params.framebuffer.w);
  glClearColor(params.background.x, params.background.y, params.background.z,
      params.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glUseProgram(glimage->program_id);
  glActiveTexture(GL_TEXTURE0 + 0);
  glBindTexture(GL_TEXTURE_2D, glimage->texture_id);
  glUniform1i(glGetUniformLocation(glimage->program_id, "txt"), 0);
  glUniform2f(glGetUniformLocation(glimage->program_id, "window_size"),
      (float)params.window.x, (float)params.window.y);
  glUniform2f(glGetUniformLocation(glimage->program_id, "image_size"),
      (float)glimage->texture_size.x, (float)glimage->texture_size.y);
  glUniform2f(glGetUniformLocation(glimage->program_id, "image_center"),
      params.center.x, params.center.y);
  glUniform1f(
      glGetUniformLocation(glimage->program_id, "image_scale"), params.scale);
  glBindBuffer(GL_ARRAY_BUFFER, glimage->texcoords_id);
  glEnableVertexAttribArray(
      glGetAttribLocation(glimage->program_id, "texcoord"));
  glVertexAttribPointer(glGetAttribLocation(glimage->program_id, "texcoord"), 2,
      GL_FLOAT, false, 0, nullptr);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glimage->triangles_id);
  glDrawElements(GL_TRIANGLES, 2 * 3, GL_UNSIGNED_INT, nullptr);
  glUseProgram(0);
  assert(glGetError() == GL_NO_ERROR);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

static const char* glscene_vertex =
    R"(
#version 330

layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
layout(location = 3) in vec4 vert_color;          // vertex color
layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

uniform mat4 shape_xform;                    // shape transform
uniform mat4 shape_xform_invtranspose;       // shape transform
uniform float shape_normal_offset;           // shape normal offset

uniform mat4 cam_xform;          // camera xform
uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
uniform mat4 cam_proj;           // camera projection

out vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
out vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
out vec4 color;                 // [to fragment shader] vertex color
out vec4 tangsp;                // [to fragment shader] vertex tangent space

// main function
void main() {
    // copy values
    pos = vert_pos;
    norm = vert_norm;
    tangsp = vert_tangsp;

    // normal offset
    if(shape_normal_offset != 0) {
        pos += shape_normal_offset * norm;
    }

    // world projection
    pos = (shape_xform * vec4(pos,1)).xyz;
    norm = (shape_xform_invtranspose * vec4(norm,0)).xyz;
    tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

    // copy other vertex properties
    texcoord = vert_texcoord;
    color = vert_color;

    // clip
    gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
}
)";

static const char* glscene_fragment =
    R"(
#version 330

float pif = 3.14159265;

uniform bool eyelight;         // eyelight shading
uniform vec3 lamb;             // ambient light
uniform int lnum;              // number of lights
uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
uniform vec3 lpos[16];         // light positions
uniform vec3 lke[16];          // light intensities

void evaluate_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
    cl = vec3(0,0,0);
    wi = vec3(0,0,0);
    if(ltype[lid] == 0) {
        // compute point light color at pos
        cl = lke[lid] / pow(length(lpos[lid]-pos),2);
        // compute light direction at pos
        wi = normalize(lpos[lid]-pos);
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

uniform int elem_type;
uniform bool elem_faceted;
uniform vec4 highlight;   // highlighted color

uniform int mat_type;          // material type
uniform vec3 mat_ke;           // material ke
uniform vec3 mat_kd;           // material kd
uniform vec3 mat_ks;           // material ks
uniform float mat_rs;          // material rs
uniform float mat_op;          // material op

uniform bool mat_ke_tex_on;    // material ke texture on
uniform sampler2D mat_ke_tex;  // material ke texture
uniform bool mat_kd_tex_on;    // material kd texture on
uniform sampler2D mat_kd_tex;  // material kd texture
uniform bool mat_ks_tex_on;    // material ks texture on
uniform sampler2D mat_ks_tex;  // material ks texture
uniform bool mat_rs_tex_on;    // material rs texture on
uniform sampler2D mat_rs_tex;  // material rs texture
uniform bool mat_op_tex_on;    // material op texture on
uniform sampler2D mat_op_tex;  // material op texture

uniform bool mat_norm_tex_on;    // material norm texture on
uniform sampler2D mat_norm_tex;  // material norm texture

uniform bool mat_double_sided;   // double sided rendering

uniform mat4 shape_xform;              // shape transform
uniform mat4 shape_xform_invtranspose; // shape transform

bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                    out vec3 kd, out vec3 ks, out float rs, out float op) {
    if(mat_type == 0) {
        ke = mat_ke;
        kd = vec3(0,0,0);
        ks = vec3(0,0,0);
        op = 1;
        return false;
    }

    ke = color.xyz * mat_ke;
    kd = color.xyz * mat_kd;
    ks = color.xyz * mat_ks;
    rs = mat_rs;
    op = color.w * mat_op;

    vec4 ke_tex = (mat_ke_tex_on) ? texture(mat_ke_tex,texcoord) : vec4(1,1,1,1);
    vec4 kd_tex = (mat_kd_tex_on) ? texture(mat_kd_tex,texcoord) : vec4(1,1,1,1);
    vec4 ks_tex = (mat_ks_tex_on) ? texture(mat_ks_tex,texcoord) : vec4(1,1,1,1);
    vec4 rs_tex = (mat_rs_tex_on) ? texture(mat_rs_tex,texcoord) : vec4(1,1,1,1);
    vec4 op_tex = (mat_op_tex_on) ? texture(mat_op_tex,texcoord) : vec4(1,1,1,1);

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

vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
    if(!mat_norm_tex_on) return norm;
    vec3 tangu = normalize((shape_xform * vec4(normalize(tangsp.xyz),0)).xyz);
    vec3 tangv = normalize(cross(norm, tangu));
    if(tangsp.w < 0) tangv = -tangv;
    vec3 texture = 2 * pow(texture(mat_norm_tex,texcoord).xyz, vec3(1/2.2)) - 1;
    // texture.y = -texture.y;
    return normalize( tangu * texture.x + tangv * texture.y + norm * texture.z );
}

in vec3 pos;                   // [from vertex shader] position in world space
in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
in vec2 texcoord;              // [from vertex shader] texcoord
in vec4 color;                 // [from vertex shader] color
in vec4 tangsp;                // [from vertex shader] tangent space

uniform vec3 cam_pos;          // camera position
uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
uniform mat4 cam_proj;           // camera projection

uniform float exposure; 
uniform float gamma;

out vec4 frag_color;      

vec3 triangle_normal(vec3 pos) {
    vec3 fdx = dFdx(pos); 
    vec3 fdy = dFdy(pos); 
    return normalize((shape_xform * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

// main
void main() {
    // view vector
    vec3 wo = normalize(cam_pos - pos);

    // prepare normals
    vec3 n;
    if(elem_faceted) {
        n = triangle_normal(pos);
    } else {
        n = normalize(norm);
    }

    // apply normal map
    n = apply_normal_map(texcoord, n, tangsp);

    // use faceforward to ensure the normals points toward us
    if(mat_double_sided) n = faceforward(n,-wo,n);

    // get material color from textures
    vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
    bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

    // exit if needed
    if(brdf_op < 0.005) discard;

    // check const color
    if(elem_type == 0) {
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
            c += pif * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
        } else {
            // accumulate ambient
            c += lamb * brdf_kd;
            // foreach light
            for(int lid = 0; lid < lnum; lid ++) {
                vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                evaluate_light(lid, pos, cl, wi);
                c += cl * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
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

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

opengl_texture::~opengl_texture() {
  if (texture_id) glDeleteTextures(1, &texture_id);
}

opengl_shape::~opengl_shape() {
  if (positions_id) glDeleteBuffers(1, &positions_id);
  if (normals_id) glDeleteBuffers(1, &normals_id);
  if (texcoords_id) glDeleteBuffers(1, &texcoords_id);
  if (colors_id) glDeleteBuffers(1, &colors_id);
  if (tangents_id) glDeleteBuffers(1, &tangents_id);
  if (points_id) glDeleteBuffers(1, &points_id);
  if (lines_id) glDeleteBuffers(1, &lines_id);
  if (triangles_id) glDeleteBuffers(1, &triangles_id);
  if (quads_id) glDeleteBuffers(1, &quads_id);
  if (edges_id) glDeleteBuffers(1, &edges_id);
}

opengl_scene::~opengl_scene() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto instance : instances) delete instance;
  for (auto texture : textures) delete texture;
  for (auto light : lights) delete light;
  if (program_id) glDeleteProgram(program_id);
  if (vertex_id) glDeleteShader(vertex_id);
  if (fragment_id) glDeleteShader(fragment_id);
  if (array_id) glDeleteVertexArrays(1, &array_id);
}

// Initialize an OpenGL scene
unique_ptr<opengl_scene> make_glscene() {
  auto glscene = make_unique<opengl_scene>();
  init_glprogram(glscene->program_id, glscene->vertex_id, glscene->fragment_id,
      glscene->array_id, glscene_vertex, glscene_fragment);
  return glscene;
}
void init_glscene(opengl_scene* glscene) {
  if (glscene->program_id) return;
  init_glprogram(glscene->program_id, glscene->vertex_id, glscene->fragment_id,
      glscene->array_id, glscene_vertex, glscene_fragment);
}
bool is_initialized(opengl_scene* glscene) { return (bool)glscene->program_id; }

// add camera
opengl_camera* add_camera(opengl_scene* scene) {
  return scene->cameras.emplace_back(new opengl_camera{});
}
void set_frame(opengl_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(opengl_camera* camera, float lens, float aspect, float film) {
  camera->lens   = lens;
  camera->aspect = aspect;
  camera->film   = film;
}
void set_nearfar(opengl_camera* camera, float near, float far) {
  camera->near = near;
  camera->far  = far;
}

// add texture
opengl_texture* add_texture(opengl_scene* scene) {
  return scene->textures.emplace_back(new opengl_texture{});
}

void set_texture(
    opengl_texture* texture, const image<vec4b>& img, bool as_srgb) {
  assert(glGetError() == GL_NO_ERROR);
  if (img.empty()) {
    glDeleteTextures(1, &texture->texture_id);
    texture->texture_id = 0;
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != img.size() || texture->nchan != 4 ||
      texture->is_srgb != as_srgb || texture->is_float == true) {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0, as_srgb ? GL_SRGB_ALPHA : GL_RGBA,
        img.size().x, img.size().y, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
        GL_UNSIGNED_BYTE, img.data());
    glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size     = img.size();
  texture->nchan    = 4;
  texture->is_srgb  = as_srgb;
  texture->is_float = false;
  assert(glGetError() == GL_NO_ERROR);
}
void set_texture(
    opengl_texture* texture, const image<vec4f>& img, bool as_float) {
  assert(glGetError() == GL_NO_ERROR);
  if (img.empty()) {
    glDeleteTextures(1, &texture->texture_id);
    texture->texture_id = 0;
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != img.size() || texture->nchan != 4 ||
      texture->is_float != as_float || texture->is_srgb == true) {
    glGenTextures(1, &texture->texture_id);
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0, as_float ? GL_RGBA32F : GL_RGBA,
        img.size().x, img.size().y, 0, GL_RGBA, GL_FLOAT, img.data());
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
        GL_FLOAT, img.data());
    glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size     = img.size();
  texture->nchan    = 4;
  texture->is_srgb  = false;
  texture->is_float = as_float;
  assert(glGetError() == GL_NO_ERROR);
}

void set_texture(
    opengl_texture* texture, const image<vec3b>& img, bool as_srgb) {
  assert(glGetError() == GL_NO_ERROR);
  if (img.empty()) {
    glDeleteTextures(1, &texture->texture_id);
    texture->texture_id = 0;
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != img.size() || texture->nchan != 3 ||
      texture->is_srgb != as_srgb || texture->is_float == true) {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0, as_srgb ? GL_SRGB : GL_RGB, img.size().x,
        img.size().y, 0, GL_RGB, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGB,
        GL_UNSIGNED_BYTE, img.data());
    glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size     = img.size();
  texture->nchan    = 3;
  texture->is_srgb  = as_srgb;
  texture->is_float = false;
  assert(glGetError() == GL_NO_ERROR);
}
void set_texture(
    opengl_texture* texture, const image<vec3f>& img, bool as_float) {
  assert(glGetError() == GL_NO_ERROR);
  if (img.empty()) {
    glDeleteTextures(1, &texture->texture_id);
    texture->texture_id = 0;
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != img.size() || texture->nchan != 3 ||
      texture->is_float != as_float || texture->is_srgb == true) {
    glGenTextures(1, &texture->texture_id);
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0, as_float ? GL_RGB32F : GL_RGB, img.size().x,
        img.size().y, 0, GL_RGB, GL_FLOAT, img.data());
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGB,
        GL_FLOAT, img.data());
    glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size     = img.size();
  texture->nchan    = 3;
  texture->is_srgb  = false;
  texture->is_float = as_float;
  assert(glGetError() == GL_NO_ERROR);
}

// add shape
opengl_shape* add_shape(opengl_scene* scene) {
  return scene->shapes.emplace_back(new opengl_shape{});
}

static void set_glshape_buffer(uint& array_id, int& array_num, bool element,
    int size, int count, const float* values) {
  if (!size || !values) {
    if (array_id) glDeleteBuffers(1, &array_id);
    array_id  = 0;
    array_num = 0;
    return;
  }
  if (!array_id) {
    init_glbuffer(array_id, element, size, count, values);
    array_num = size;
  } else if (array_num != size) {
    glDeleteBuffers(1, &array_id);
    init_glbuffer(array_id, element, size, count, values);
    array_num = size;
  } else {
    update_glbuffer(array_id, element, size, 3, values);
  }
}
static void set_glshape_buffer(uint& array_id, int& array_num, bool element,
    int size, int count, const int* values) {
  if (!size || !values) {
    if (array_id) glDeleteBuffers(1, &array_id);
    array_id  = 0;
    array_num = 0;
    return;
  }
  if (!array_id) {
    init_glbuffer(array_id, element, size, count, values);
    array_num = size;
  } else if (array_num != size) {
    glDeleteBuffers(1, &array_id);
    init_glbuffer(array_id, element, size, count, values);
    array_num = size;
  } else {
    update_glbuffer(array_id, element, size, 3, values);
  }
}

void set_points(opengl_shape* shape, const vector<int>& points) {
  set_glshape_buffer(shape->points_id, shape->points_num, true, points.size(),
      1, (const int*)points.data());
}
void set_lines(opengl_shape* shape, const vector<vec2i>& lines) {
  set_glshape_buffer(shape->lines_id, shape->lines_num, true, lines.size(), 2,
      (const int*)lines.data());
}
void set_triangles(opengl_shape* shape, const vector<vec3i>& triangles) {
  set_glshape_buffer(shape->triangles_id, shape->triangles_num, true,
      triangles.size(), 3, (const int*)triangles.data());
}
void set_quads(opengl_shape* shape, const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  set_glshape_buffer(shape->quads_id, shape->quads_num, true, triangles.size(),
      3, (const int*)triangles.data());
}
void set_positions(opengl_shape* shape, const vector<vec3f>& positions) {
  set_glshape_buffer(shape->positions_id, shape->positions_num, false,
      positions.size(), 3, (const float*)positions.data());
}
void set_normals(opengl_shape* shape, const vector<vec3f>& normals) {
  set_glshape_buffer(shape->normals_id, shape->normals_num, false,
      normals.size(), 3, (const float*)normals.data());
}
void set_texcoords(opengl_shape* shape, const vector<vec2f>& texcoords) {
  set_glshape_buffer(shape->texcoords_id, shape->texcoords_num, false,
      texcoords.size(), 2, (const float*)texcoords.data());
}
void set_colors(opengl_shape* shape, const vector<vec3f>& colors) {
  set_glshape_buffer(shape->colors_id, shape->colors_num, false, colors.size(),
      3, (const float*)colors.data());
}
void set_tangents(opengl_shape* shape, const vector<vec4f>& tangents) {
  set_glshape_buffer(shape->tangents_id, shape->tangents_num, false,
      tangents.size(), 4, (const float*)tangents.data());
}

// add object
opengl_object* add_object(opengl_scene* scene) {
  return scene->objects.emplace_back(new opengl_object{});
}
void set_frame(opengl_object* object, const frame3f& frame) {
  object->frame = frame;
}
void set_shape(opengl_object* object, opengl_shape* shape) {
  object->shape = shape;
}
void set_material(opengl_object* object, opengl_material* material) {
  object->material = material;
}
void set_instance(opengl_object* object, opengl_instance* instance) {
  object->instance = instance;
}
void set_hidden(opengl_object* object, bool hidden) { object->hidden = hidden; }
void set_highlighted(opengl_object* object, bool highlighted) {
  object->highlighted = highlighted;
}

// add instance
opengl_instance* add_instance(opengl_scene* scene) {
  return scene->instances.emplace_back(new opengl_instance{});
}
void set_frames(opengl_instance* instance, const vector<frame3f>& frames) {
  // TODO: instances
}

// add material
opengl_material* add_material(opengl_scene* scene) {
  return scene->materials.emplace_back(new opengl_material{});
}
void set_emission(opengl_material* material, const vec3f& emission,
    opengl_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    opengl_material* material, const vec3f& color, opengl_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    opengl_material* material, float specular, opengl_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_roughness(
    opengl_material* material, float roughness, opengl_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    opengl_material* material, float opacity, opengl_texture* opacity_tex) {
  material->opacity = opacity;
}
void set_metallic(
    opengl_material* material, float metallic, opengl_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_normalmap(opengl_material* material, opengl_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// add light
opengl_light* add_light(opengl_scene* scene) {
  return scene->lights.emplace_back(new opengl_light{});
}
void set_light(opengl_light* light, const vec3f& position,
    const vec3f& emission, bool directional) {
  light->position = position;
  light->emission = emission;
  light->type     = directional ? 1 : 0;
}
void clear_lights(opengl_scene* scene) {
  for (auto light : scene->lights) delete light;
  scene->lights.clear();
}
bool has_max_lights(opengl_scene* scene) { return scene->lights.size() >= 16; }

// Draw a shape
void draw_globject(opengl_scene* glscene, opengl_object* object,
    const draw_glscene_params& params) {
  if (object->hidden) return;

  auto instance_xform     = mat4f(object->frame);
  auto instance_inv_xform = transpose(
      mat4f(inverse(object->frame, params.non_rigid_frames)));
  glUniformMatrix4fv(glGetUniformLocation(glscene->program_id, "shape_xform"),
      1, false, &instance_xform.x.x);
  glUniformMatrix4fv(
      glGetUniformLocation(glscene->program_id, "shape_xform_invtranspose"), 1,
      false, &instance_inv_xform.x.x);
  glUniform1f(
      glGetUniformLocation(glscene->program_id, "shape_normal_offset"), 0.0f);
  if (object->highlighted) {
    glUniform4f(
        glGetUniformLocation(glscene->program_id, "highlight"), 1, 1, 0, 1);
  } else {
    glUniform4f(
        glGetUniformLocation(glscene->program_id, "highlight"), 0, 0, 0, 0);
  }

  auto material = object->material;
  auto mtype    = 2;
  glUniform1i(glGetUniformLocation(glscene->program_id, "mat_type"), mtype);
  glUniform3f(glGetUniformLocation(glscene->program_id, "mat_ke"),
      material->emission.x, material->emission.y, material->emission.z);
  glUniform3f(glGetUniformLocation(glscene->program_id, "mat_kd"),
      material->color.x, material->color.y, material->color.z);
  glUniform3f(glGetUniformLocation(glscene->program_id, "mat_ks"),
      material->metallic, material->metallic, material->metallic);
  glUniform1f(
      glGetUniformLocation(glscene->program_id, "mat_rs"), material->roughness);
  glUniform1f(
      glGetUniformLocation(glscene->program_id, "mat_op"), material->opacity);
  glUniform1i(glGetUniformLocation(glscene->program_id, "mat_double_sided"),
      (int)params.double_sided);
  if (material->emission_tex) {
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_2D, material->emission_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ke_tex"), 0);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ke_tex_on"), 1);
  } else {
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ke_tex_on"), 0);
  }
  if (material->color_tex) {
    glActiveTexture(GL_TEXTURE0 + 1);
    glBindTexture(GL_TEXTURE_2D, material->color_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_kd_tex"), 1);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_kd_tex_on"), 1);
  } else {
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_kd_tex_on"), 0);
  }
  if (material->metallic_tex) {
    glActiveTexture(GL_TEXTURE0 + 2);
    glBindTexture(GL_TEXTURE_2D, material->specular_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ks_tex"), 2);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ks_tex_on"), 1);
  } else {
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_ks_tex_on"), 0);
  }
  if (material->roughness_tex) {
    glActiveTexture(GL_TEXTURE0 + 3);
    glBindTexture(GL_TEXTURE_2D, material->roughness_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_rs_tex"), 3);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_rs_tex_on"), 1);
  } else {
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_rs_tex_on"), 0);
  }
  if (material->opacity_tex) {
    glActiveTexture(GL_TEXTURE0 + 3);
    glBindTexture(GL_TEXTURE_2D, material->opacity_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_op_tex"), 3);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_op_tex_on"), 1);
  } else {
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_op_tex_on"), 0);
  }
  if (material->normal_tex) {
    glActiveTexture(GL_TEXTURE0 + 4);
    glBindTexture(GL_TEXTURE_2D, material->normal_tex->texture_id);
    glUniform1i(glGetUniformLocation(glscene->program_id, "mat_norm_tex"), 4);
    glUniform1i(
        glGetUniformLocation(glscene->program_id, "mat_norm_tex_on"), 1);
  } else {
    glUniform1i(
        glGetUniformLocation(glscene->program_id, "mat_norm_tex_on"), 0);
  }

  auto shape = object->shape;
  glUniform1i(glGetUniformLocation(glscene->program_id, "elem_faceted"),
      (int)!shape->normals_id);
  if (shape->positions_id) {
    glBindBuffer(GL_ARRAY_BUFFER, shape->positions_id);
    glEnableVertexAttribArray(
        glGetAttribLocation(glscene->program_id, "vert_pos"));
    glVertexAttribPointer(glGetAttribLocation(glscene->program_id, "vert_pos"),
        3, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib3f(
        glGetAttribLocation(glscene->program_id, "vert_pos"), 0, 0, 0);
  }
  if (shape->normals_id) {
    glBindBuffer(GL_ARRAY_BUFFER, shape->normals_id);
    glEnableVertexAttribArray(
        glGetAttribLocation(glscene->program_id, "vert_norm"));
    glVertexAttribPointer(glGetAttribLocation(glscene->program_id, "vert_norm"),
        3, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib3f(
        glGetAttribLocation(glscene->program_id, "vert_norm"), 0, 0, 0);
  }
  if (shape->texcoords_id) {
    glBindBuffer(GL_ARRAY_BUFFER, shape->texcoords_id);
    glEnableVertexAttribArray(
        glGetAttribLocation(glscene->program_id, "vert_texcoord"));
    glVertexAttribPointer(
        glGetAttribLocation(glscene->program_id, "vert_texcoord"), 2, GL_FLOAT,
        false, 0, nullptr);
  } else {
    glVertexAttrib2f(
        glGetAttribLocation(glscene->program_id, "vert_texcoord"), 0, 0);
  }
  if (shape->colors_id) {
    glBindBuffer(GL_ARRAY_BUFFER, shape->colors_id);
    glEnableVertexAttribArray(
        glGetAttribLocation(glscene->program_id, "vert_color"));
    glVertexAttribPointer(
        glGetAttribLocation(glscene->program_id, "vert_color"), 4, GL_FLOAT,
        false, 0, nullptr);
  } else {
    glVertexAttrib4f(
        glGetAttribLocation(glscene->program_id, "vert_color"), 1, 1, 1, 1);
  }
  if (shape->tangents_id) {
    glBindBuffer(GL_ARRAY_BUFFER, shape->tangents_id);
    glEnableVertexAttribArray(
        glGetAttribLocation(glscene->program_id, "vert_tangsp"));
    glVertexAttribPointer(
        glGetAttribLocation(glscene->program_id, "vert_tangsp"), 4, GL_FLOAT,
        false, 0, nullptr);
  } else {
    glVertexAttrib4f(
        glGetAttribLocation(glscene->program_id, "vert_tangsp"), 0, 0, 1, 1);
  }

  if (shape->points_id) {
    glUniform1i(glGetUniformLocation(glscene->program_id, "elem_type"), 1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape->points_id);
    glDrawElements(GL_POINTS, shape->points_num, GL_UNSIGNED_INT, nullptr);
  }
  if (shape->lines_id) {
    glUniform1i(glGetUniformLocation(glscene->program_id, "elem_type"), 2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape->lines_id);
    glDrawElements(GL_LINES, shape->lines_num * 2, GL_UNSIGNED_INT, nullptr);
  }
  if (shape->triangles_id) {
    glUniform1i(glGetUniformLocation(glscene->program_id, "elem_type"), 3);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape->triangles_id);
    glDrawElements(
        GL_TRIANGLES, shape->triangles_num * 3, GL_UNSIGNED_INT, nullptr);
  }
  if (shape->quads_id) {
    glUniform1i(glGetUniformLocation(glscene->program_id, "elem_type"), 3);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape->quads_id);
    glDrawElements(
        GL_TRIANGLES, shape->quads_num * 3, GL_UNSIGNED_INT, nullptr);
  }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_glerror();
        set_gluniform(glscene->program, "mtype"), 0);
        glUniform3f(glGetUniformLocation(glscene->program, "ke"), 0, 0, 0);
        set_gluniform(glscene->program, "op"), shape->op);
        set_gluniform(glscene->program, "shp_normal_offset"), 0.01f);
        check_glerror();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_glerror();
    }
#endif
  if (params.edges) throw std::runtime_error("edges are momentarily disabled");
  // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
}

// Display a scene
void draw_glscene(opengl_scene* glscene, const vec4i& viewport,
    const draw_glscene_params& params) {
  auto& glcamera      = glscene->cameras.at(params.camera);
  auto  camera_aspect = (float)viewport.z / (float)viewport.w;
  auto  camera_yfov =
      camera_aspect >= 0
          ? (2 * atan(glcamera->film / (camera_aspect * 2 * glcamera->lens)))
          : (2 * atan(glcamera->film / (2 * glcamera->lens)));
  auto camera_view = mat4f(inverse(glcamera->frame));
  auto camera_proj = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  glClearColor(params.background.x, params.background.y, params.background.z,
      params.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);

  glUseProgram(glscene->program_id);
  glUniform3f(glGetUniformLocation(glscene->program_id, "cam_pos"),
      glcamera->frame.o.x, glcamera->frame.o.y, glcamera->frame.o.z);
  glUniformMatrix4fv(glGetUniformLocation(glscene->program_id, "cam_xform_inv"),
      1, false, &camera_view.x.x);
  glUniformMatrix4fv(glGetUniformLocation(glscene->program_id, "cam_proj"), 1,
      false, &camera_proj.x.x);
  glUniform1i(glGetUniformLocation(glscene->program_id, "eyelight"),
      (int)params.eyelight);
  glUniform1f(
      glGetUniformLocation(glscene->program_id, "exposure"), params.exposure);
  glUniform1f(glGetUniformLocation(glscene->program_id, "gamma"), params.gamma);

  if (!params.eyelight) {
    glUniform3f(glGetUniformLocation(glscene->program_id, "lamb"), 0, 0, 0);
    glUniform1i(glGetUniformLocation(glscene->program_id, "lnum"),
        (int)glscene->lights.size());
    for (auto i = 0; i < glscene->lights.size(); i++) {
      auto is = std::to_string(i);
      glUniform3f(glGetUniformLocation(
                      glscene->program_id, ("lpos[" + is + "]").c_str()),
          glscene->lights[i]->position.x, glscene->lights[i]->position.y,
          glscene->lights[i]->position.z);
      glUniform3f(glGetUniformLocation(
                      glscene->program_id, ("lke[" + is + "]").c_str()),
          glscene->lights[i]->emission.x, glscene->lights[i]->emission.y,
          glscene->lights[i]->emission.z);
      glUniform1i(glGetUniformLocation(
                      glscene->program_id, ("ltype[" + is + "]").c_str()),
          (int)glscene->lights[i]->type);
    }
  }

  if (params.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  for (auto object : glscene->objects) {
    draw_globject(glscene, object, params);
  }

  glUseProgram(0);
  if (params.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

static void draw_glwindow(opengl_window* win) {
  glClearColor(win->background.x, win->background.y, win->background.z,
      win->background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (win->draw_cb) win->draw_cb(win, win->input);
  if (win->widgets_cb) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    auto window = zero2i;
    glfwGetWindowSize(win->win, &window.x, &window.y);
    if (win->widgets_left) {
      ImGui::SetNextWindowPos({0, 0});
      ImGui::SetNextWindowSize({(float)win->widgets_width, (float)window.y});
    } else {
      ImGui::SetNextWindowPos({(float)(window.x - win->widgets_width), 0});
      ImGui::SetNextWindowSize({(float)win->widgets_width, (float)window.y});
    }
    ImGui::SetNextWindowCollapsed(false);
    ImGui::SetNextWindowBgAlpha(1);
    if (ImGui::Begin(win->title.c_str(), nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings)) {
      draw_glmessages(win);
      win->widgets_cb(win, win->input);
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  }
  glfwSwapBuffers(win->win);
}

unique_ptr<opengl_window> make_glwindow(const vec2i& size, const string& title,
    bool widgets, int widgets_width, bool widgets_left) {
  auto win = make_unique<opengl_window>();

  // init glfw
  if (!glfwInit())
    throw std::runtime_error("cannot initialize windowing system");
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // create window
  win->title = title;
  win->win = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
  if (!win->win) throw std::runtime_error("cannot initialize windowing system");
  glfwMakeContextCurrent(win->win);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowUserPointer(win->win, win.get());

  // set callbacks
  glfwSetWindowRefreshCallback(win->win, [](GLFWwindow* glfw) {
    auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
    draw_glwindow(win);
  });
  glfwSetDropCallback(
      win->win, [](GLFWwindow* glfw, int num, const char** paths) {
        auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
        if (win->drop_cb) {
          auto pathv = vector<string>();
          for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
          win->drop_cb(win, pathv, win->input);
        }
      });
  glfwSetKeyCallback(win->win,
      [](GLFWwindow* glfw, int key, int scancode, int action, int mods) {
        auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
        if (win->key_cb) win->key_cb(win, key, (bool)action, win->input);
      });
  glfwSetMouseButtonCallback(
      win->win, [](GLFWwindow* glfw, int button, int action, int mods) {
        auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
        if (win->click_cb)
          win->click_cb(
              win, button == GLFW_MOUSE_BUTTON_LEFT, (bool)action, win->input);
      });
  glfwSetScrollCallback(
      win->win, [](GLFWwindow* glfw, double xoffset, double yoffset) {
        auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
        if (win->scroll_cb) win->scroll_cb(win, (float)yoffset, win->input);
      });
  glfwSetWindowSizeCallback(
      win->win, [](GLFWwindow* glfw, int width, int height) {
        auto win = (opengl_window*)glfwGetWindowUserPointer(glfw);
        glfwGetWindowSize(
            win->win, &win->input.window_size.x, &win->input.window_size.y);
        if (win->widgets_width) win->input.window_size.x -= win->widgets_width;
        glfwGetFramebufferSize(win->win, &win->input.framebuffer_viewport.z,
            &win->input.framebuffer_viewport.w);
        win->input.framebuffer_viewport.x = 0;
        win->input.framebuffer_viewport.y = 0;
        if (win->widgets_width) {
          auto win_size = zero2i;
          glfwGetWindowSize(win->win, &win_size.x, &win_size.y);
          auto offset = (int)(win->widgets_width *
                              (float)win->input.framebuffer_viewport.z /
                              win_size.x);
          win->input.framebuffer_viewport.z -= offset;
          if (win->widgets_left) win->input.framebuffer_viewport.x += offset;
        }
      });

  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error("cannot initialize OpenGL extensions");

  // widgets
  if (widgets) {
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename       = nullptr;
    ImGui::GetStyle().WindowRounding = 0;
    ImGui_ImplGlfw_InitForOpenGL(win->win, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    ImGui::StyleColorsDark();
    win->widgets_width = widgets_width;
    win->widgets_left  = widgets_left;
  }

  return win;
}

void clear_glwindow(opengl_window* win) {
  glfwDestroyWindow(win->win);
  glfwTerminate();
  win->win = nullptr;
}

// Run loop
void run_ui(opengl_window* win) {
  while (!glfwWindowShouldClose(win->win)) {
    // update input
    win->input.mouse_last = win->input.mouse_pos;
    auto mouse_posx = 0.0, mouse_posy = 0.0;
    glfwGetCursorPos(win->win, &mouse_posx, &mouse_posy);
    win->input.mouse_pos = vec2f{(float)mouse_posx, (float)mouse_posy};
    if (win->widgets_width && win->widgets_left)
      win->input.mouse_pos.x -= win->widgets_width;
    win->input.mouse_left = glfwGetMouseButton(
                                win->win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    win->input.mouse_right =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    win->input.modifier_alt =
        glfwGetKey(win->win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
    win->input.modifier_shift =
        glfwGetKey(win->win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    win->input.modifier_ctrl =
        glfwGetKey(win->win, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
    glfwGetWindowSize(
        win->win, &win->input.window_size.x, &win->input.window_size.y);
    if (win->widgets_width) win->input.window_size.x -= win->widgets_width;
    glfwGetFramebufferSize(win->win, &win->input.framebuffer_viewport.z,
        &win->input.framebuffer_viewport.w);
    win->input.framebuffer_viewport.x = 0;
    win->input.framebuffer_viewport.y = 0;
    if (win->widgets_width) {
      auto win_size = zero2i;
      glfwGetWindowSize(win->win, &win_size.x, &win_size.y);
      auto offset = (int)(win->widgets_width *
                          (float)win->input.framebuffer_viewport.z /
                          win_size.x);
      win->input.framebuffer_viewport.z -= offset;
      if (win->widgets_left) win->input.framebuffer_viewport.x += offset;
    }
    if (win->widgets_width) {
      auto io                   = &ImGui::GetIO();
      win->input.widgets_active = io->WantTextInput || io->WantCaptureMouse ||
                                  io->WantCaptureKeyboard;
    }

    // time
    win->input.clock_last = win->input.clock_now;
    win->input.clock_now =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    win->input.time_now = (double)win->input.clock_now / 1000000000.0;
    win->input.time_delta =
        (double)(win->input.clock_now - win->input.clock_last) / 1000000000.0;

    // update ui
    if (win->uiupdate_cb) win->uiupdate_cb(win, win->input);

    // update
    if (win->update_cb) win->update_cb(win, win->input);

    // draw
    draw_glwindow(win);

    // event hadling
    glfwPollEvents();
  }
}

void set_draw_glcallback(opengl_window* win, draw_glcallback cb) {
  win->draw_cb = cb;
}
void set_widgets_glcallback(opengl_window* win, widgets_glcallback cb) {
  win->widgets_cb = cb;
}
void set_drop_glcallback(opengl_window* win, drop_glcallback drop_cb) {
  win->drop_cb = drop_cb;
}
void set_key_glcallback(opengl_window* win, key_glcallback cb) {
  win->key_cb = cb;
}
void set_click_glcallback(opengl_window* win, click_glcallback cb) {
  win->click_cb = cb;
}
void set_scroll_glcallback(opengl_window* win, scroll_glcallback cb) {
  win->scroll_cb = cb;
}
void set_uiupdate_glcallback(opengl_window* win, uiupdate_glcallback cb) {
  win->uiupdate_cb = cb;
}
void set_update_glcallback(opengl_window* win, update_glcallback cb) {
  win->update_cb = cb;
}

void set_close(opengl_window* win, bool close) {
  glfwSetWindowShouldClose(win->win, close ? GLFW_TRUE : GLFW_FALSE);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

void init_glwidgets(opengl_window* win, int width, bool left) {
  // init widgets
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(win->win, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();
  win->widgets_width = width;
  win->widgets_left  = left;
}

bool begin_glheader(opengl_window* win, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_glheader(opengl_window* win) { ImGui::PopID(); }

void open_glmodal(opengl_window* win, const char* lbl) {
  ImGui::OpenPopup(lbl);
}
void clear_glmodal(opengl_window* win) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(opengl_window* win, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(opengl_window* win) { ImGui::EndPopup(); }
bool is_glmodal_open(opengl_window* win, const char* lbl) {
  return ImGui::IsPopupOpen(lbl);
}

bool draw_glmessage(
    opengl_window* win, const char* lbl, const string& message) {
  if (ImGui::BeginPopupModal(lbl)) {
    auto open = true;
    ImGui::Text("%s", message.c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      open = false;
    }
    ImGui::EndPopup();
    return open;
  } else {
    return false;
  }
}

std::deque<string> _message_queue = {};
std::mutex         _message_mutex;
void               push_glmessage(opengl_window* win, const string& message) {
  std::lock_guard lock(_message_mutex);
  _message_queue.push_back(message);
}
bool draw_glmessages(opengl_window* win) {
  std::lock_guard lock(_message_mutex);
  if (_message_queue.empty()) return false;
  if (!is_glmodal_open(win, "<message>")) {
    open_glmodal(win, "<message>");
    return true;
  } else if (ImGui::BeginPopupModal("<message>")) {
    ImGui::Text("%s", _message_queue.front().c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      _message_queue.pop_front();
    }
    ImGui::EndPopup();
    return true;
  } else {
    return false;
  }
}

// Utility to normalize a path
static inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)

    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get extension (not including '.').
static string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

struct filedialog_state {
  string                     dirname       = "";
  string                     filename      = "";
  vector<pair<string, bool>> entries       = {};
  bool                       save          = false;
  bool                       remove_hidden = true;
  string                     filter        = "";
  vector<string>             extensions    = {};

  filedialog_state() {}
  filedialog_state(const string& dirname, const string& filename, bool save,
      const string& filter) {
    this->save = save;
    set_filter(filter);
    set_dirname(dirname);
    set_filename(filename);
  }
  void set_dirname(const string& name) {
    dirname = name;
    dirname = normalize_path(dirname);
    if (dirname == "") dirname = "./";
    if (dirname.back() != '/') dirname += '/';
    refresh();
  }
  void set_filename(const string& name) {
    filename = name;
    check_filename();
  }
  void set_filter(const string& flt) {
    auto globs = vector<string>{""};
    for (auto i = 0; i < flt.size(); i++) {
      if (flt[i] == ';') {
        globs.push_back("");
      } else {
        globs.back() += flt[i];
      }
    }
    filter = "";
    extensions.clear();
    for (auto pattern : globs) {
      if (pattern == "") continue;
      auto ext = get_extension(pattern);
      if (ext != "") {
        extensions.push_back(ext);
        filter += (filter == "") ? ("*." + ext) : (";*." + ext);
      }
    }
  }
  void check_filename() {
    if (filename.empty()) return;
    auto ext = get_extension(filename);
    if (std::find(extensions.begin(), extensions.end(), ext) ==
        extensions.end()) {
      filename = "";
      return;
    }
    if (!save && !exists_file(dirname + filename)) {
      filename = "";
      return;
    }
  }
  void select_entry(int idx) {
    if (entries[idx].second) {
      set_dirname(dirname + entries[idx].first);
    } else {
      set_filename(entries[idx].first);
    }
  }

  void refresh() {
    entries.clear();
    cf_dir_t dir;
    cf_dir_open(&dir, dirname.c_str());
    while (dir.has_next) {
      cf_file_t file;
      cf_read_file(&dir, &file);
      cf_dir_next(&dir);
      if (remove_hidden && file.name[0] == '.') continue;
      if (file.is_dir) {
        entries.push_back({file.name + "/"s, true});
      } else {
        entries.push_back({file.name, false});
      }
    }
    cf_dir_close(&dir);
    std::sort(entries.begin(), entries.end(), [](auto& a, auto& b) {
      if (a.second == b.second) return a.first < b.first;
      return a.second;
    });
  }

  string get_path() const { return dirname + filename; }
  bool   exists_file(const string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
  }
};
bool draw_glfiledialog(opengl_window* win, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter) {
  static auto states = unordered_map<string, filedialog_state>{};
  ImGui::SetNextWindowSize({500, 300}, ImGuiCond_FirstUseEver);
  if (ImGui::BeginPopupModal(lbl)) {
    if (states.find(lbl) == states.end()) {
      states[lbl] = filedialog_state{dirname, filename, save, filter};
    }
    auto& state = states.at(lbl);
    char  dir_buffer[1024];
    strcpy(dir_buffer, state.dirname.c_str());
    if (ImGui::InputText("dir", dir_buffer, sizeof(dir_buffer))) {
      state.set_dirname(dir_buffer);
    }
    auto current_item = -1;
    if (ImGui::ListBox(
            "entries", &current_item,
            [](void* data, int idx, const char** out_text) -> bool {
              auto& state = *(filedialog_state*)data;
              *out_text   = state.entries[idx].first.c_str();
              return true;
            },
            &state, (int)state.entries.size())) {
      state.select_entry(current_item);
    }
    char file_buffer[1024];
    strcpy(file_buffer, state.filename.c_str());
    if (ImGui::InputText("file", file_buffer, sizeof(file_buffer))) {
      state.set_filename(file_buffer);
    }
    char filter_buffer[1024];
    strcpy(filter_buffer, state.filter.c_str());
    if (ImGui::InputText("filter", filter_buffer, sizeof(filter_buffer))) {
      state.set_filter(filter_buffer);
    }
    auto ok = false, exit = false;
    if (ImGui::Button("Ok")) {
      path = state.dirname + state.filename;
      ok   = true;
      exit = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      exit = true;
    }
    if (exit) {
      ImGui::CloseCurrentPopup();
      states.erase(lbl);
    }
    ImGui::EndPopup();
    return ok;
  } else {
    return false;
  }
}
bool draw_glfiledialog_button(opengl_window* win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter) {
  if (is_glmodal_open(win, lbl)) {
    return draw_glfiledialog(win, lbl, path, save, dirname, filename, filter);
  } else {
    if (draw_glbutton(win, button_lbl, button_active)) {
      open_glmodal(win, lbl);
    }
    return false;
  }
}

bool draw_glbutton(opengl_window* win, const char* lbl, bool enabled) {
  if (enabled) {
    return ImGui::Button(lbl);
  } else {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    auto ok = ImGui::Button(lbl);
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
    return ok;
  }
}

void draw_gllabel(opengl_window* win, const char* lbl, const string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}

void draw_glseparator(opengl_window* win) { ImGui::Separator(); }

void continue_glline(opengl_window* win) { ImGui::SameLine(); }

bool draw_gltextinput(opengl_window* win, const char* lbl, string& value) {
  char buffer[4096];
  auto num = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer, sizeof(buffer));
  if (edited) value = buffer;
  return edited;
}

bool draw_glslider(
    opengl_window* win, const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_glslider(
    opengl_window* win, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_glslider(
    opengl_window* win, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_gldragger(opengl_window* win, const char* lbl, float& value,
    float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec2f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec3f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec4f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_gldragger(opengl_window* win, const char* lbl, int& value,
    float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec2i& value,
    float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec3i& value,
    float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(opengl_window* win, const char* lbl, vec4i& value,
    float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_glcheckbox(opengl_window* win, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}

bool draw_glcoloredit(opengl_window* win, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_glcoloredit(opengl_window* win, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_glhdrcoloredit(opengl_window* win, const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_glhdrcoloredit(opengl_window* win, const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(xyz(color));
  if (scale > 1) {
    xyz(color) /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    xyz(value) = xyz(color) * exp2(exposure);
    value.w    = color.w;
    return true;
  } else {
    return false;
  }
}

bool draw_glcombobox(opengl_window* win, const char* lbl, int& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, labels[value].c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(opengl_window* win, const char* lbl, string& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(opengl_window* win, const char* lbl, int& idx, int num,
    const std::function<const char*(int)>& labels, bool include_null) {
  if (num <= 0) idx = -1;
  if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx) : "<none>")) return false;
  auto old_idx = idx;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", idx < 0)) idx = -1;
    if (idx < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < num; i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels(i), idx == i)) idx = i;
    if (idx == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return idx != old_idx;
}

void draw_glprogressbar(opengl_window* win, const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_glhistogram(
    opengl_window* win, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_glhistogram(
    opengl_window* win, const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_glhistogram(
    opengl_window* win, const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_glhistogram(
    opengl_window* win, const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_glhistogram(
    opengl_window* win, const char* lbl, const vector<vec4f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " w"s).c_str(), (const float*)values.data() + 3,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
}

// https://github.com/ocornut/imgui/issues/300
struct ImGuiAppLog {
  ImGuiTextBuffer Buf;
  ImGuiTextFilter Filter;
  ImVector<int>   LineOffsets;  // Index to lines offset
  bool            ScrollToBottom;

  void Clear() {
    Buf.clear();
    LineOffsets.clear();
  }

  void AddLog(const char* msg, const char* lbl) {
    int old_size = Buf.size();
    Buf.appendf("[%s] %s\n", lbl, msg);
    for (int new_size = Buf.size(); old_size < new_size; old_size++)
      if (Buf[old_size] == '\n') LineOffsets.push_back(old_size);
    ScrollToBottom = true;
  }

  void Draw() {
    if (ImGui::Button("Clear")) Clear();
    ImGui::SameLine();
    bool copy = ImGui::Button("Copy");
    ImGui::SameLine();
    Filter.Draw("Filter", -100.0f);
    ImGui::Separator();
    ImGui::BeginChild("scrolling");
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
    if (copy) ImGui::LogToClipboard();

    if (Filter.IsActive()) {
      const char* buf_begin = Buf.begin();
      const char* line      = buf_begin;
      for (int line_no = 0; line != NULL; line_no++) {
        const char* line_end = (line_no < LineOffsets.Size)
                                   ? buf_begin + LineOffsets[line_no]
                                   : NULL;
        if (Filter.PassFilter(line, line_end))
          ImGui::TextUnformatted(line, line_end);
        line = line_end && line_end[1] ? line_end + 1 : NULL;
      }
    } else {
      ImGui::TextUnformatted(Buf.begin());
    }

    if (ScrollToBottom) ImGui::SetScrollHere(1.0f);
    ScrollToBottom = false;
    ImGui::PopStyleVar();
    ImGui::EndChild();
  }
  void Draw(const char* title, bool* p_opened = NULL) {
    ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(title, p_opened);
    Draw();
    ImGui::End();
  }
};

std::mutex  _log_mutex;
ImGuiAppLog _log_widget;
void        log_glinfo(opengl_window* win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_glerror(opengl_window* win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_gllogs(opengl_window* win) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void draw_gllog(opengl_window* win) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

}  // namespace yocto
