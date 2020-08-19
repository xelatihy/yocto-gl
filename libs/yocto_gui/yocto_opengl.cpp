//
// Utilities to use OpenGL 3, GLFW and ImGui.
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

#include "yocto_opengl.h"

#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "ext/glad/glad.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unordered_map;
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL OPENGL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

bool init_ogl(string& error) {
  if (!gladLoadGL()) {
    error = "Cannot initialize OpenGL context.";
    return false;
  }
  return true;
}

void assert_ogl_error() { assert(glGetError() == GL_NO_ERROR); }

bool check_ogl_error(string& error) {
  if (glGetError() != GL_NO_ERROR) {
    error = "";
    return false;
  }
  return true;
}

void clear_ogl_framebuffer(const vec4f& color, bool clear_depth) {
  glClearColor(color.x, color.y, color.z, color.w);
  if (clear_depth) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
  } else {
    glClear(GL_COLOR_BUFFER_BIT);
  }
}

void set_ogl_viewport(const vec4i& viewport) {
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);
}

void set_ogl_viewport(const vec2i& viewport) {
  glViewport(0, 0, viewport.x, viewport.y);
}

void set_ogl_wireframe(bool enabled) {
  if (enabled)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void set_ogl_blending(bool enabled) {
  if (enabled) {
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
  } else {
    glDisable(GL_BLEND);
  }
}

void set_ogl_point_size(int size) { glPointSize(size); }

void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const byte* img, bool as_srgb, bool linear, bool mipmap) {
  static auto sformat = unordered_map<int, uint>{
      {1, GL_SRGB},
      {2, GL_SRGB},
      {3, GL_SRGB},
      {4, GL_SRGB_ALPHA},
  };
  static auto iformat = unordered_map<int, uint>{
      {1, GL_RGB},
      {2, GL_RGB},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  static auto cformat = unordered_map<int, uint>{
      {1, GL_RED},
      {2, GL_RG},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  assert_ogl_error();
  if (size == zero2i) {
    clear_texture(texture);
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != size || texture->nchannels != nchannels ||
      texture->is_srgb != as_srgb || texture->is_float == true ||
      texture->linear != linear || texture->mipmap != mipmap) {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0,
        as_srgb ? sformat.at(nchannels) : iformat.at(nchannels), size.x, size.y,
        0, cformat.at(nchannels), GL_UNSIGNED_BYTE, img);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
        cformat.at(nchannels), GL_UNSIGNED_BYTE, img);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size      = size;
  texture->nchannels = nchannels;
  texture->is_srgb   = as_srgb;
  texture->is_float  = false;
  texture->linear    = linear;
  texture->mipmap    = mipmap;
  assert_ogl_error();
}

void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const float* img, bool as_float, bool linear, bool mipmap) {
  static auto fformat = unordered_map<int, uint>{
      {1, GL_RGB16F},
      {2, GL_RGB16F},
      {3, GL_RGB16F},
      {4, GL_RGBA32F},
  };
  static auto iformat = unordered_map<int, uint>{
      {1, GL_RGB},
      {2, GL_RGB},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  static auto cformat = unordered_map<int, uint>{
      {1, GL_RED},
      {2, GL_RG},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  assert_ogl_error();
  if (size == zero2i) {
    clear_texture(texture);
    return;
  }
  if (!texture->texture_id) glGenTextures(1, &texture->texture_id);
  if (texture->size != size || texture->nchannels != nchannels ||
      texture->is_float != as_float || texture->is_srgb == true ||
      texture->linear != linear || texture->mipmap != mipmap) {
    glGenTextures(1, &texture->texture_id);
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0,
        as_float ? fformat.at(nchannels) : iformat.at(nchannels), size.x,
        size.y, 0, iformat.at(nchannels), GL_FLOAT, img);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
        iformat.at(nchannels), GL_FLOAT, img);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  texture->size      = size;
  texture->nchannels = nchannels;
  texture->is_srgb   = false;
  texture->is_float  = as_float;
  texture->linear    = linear;
  texture->mipmap    = mipmap;
  assert_ogl_error();
}

// check if texture is initialized
bool is_initialized(ogl_texture* texture) { return texture->texture_id != 0; }

// clear texture
void clear_texture(ogl_texture* texture) {
  if (texture->texture_id) glDeleteTextures(1, &texture->texture_id);
  texture->texture_id = 0;
  texture->size       = {0, 0};
  texture->nchannels  = 0;
  texture->is_srgb    = false;
  texture->is_float   = false;
  texture->linear     = false;
  texture->mipmap     = false;
}

void set_texture(ogl_texture* texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 4, (const byte*)img.data(), as_srgb,
      linear, mipmap);
}
void set_texture(ogl_texture* texture, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 4, (const float*)img.data(), as_float,
      linear, mipmap);
}

void set_texture(ogl_texture* texture, const image<vec3b>& img, bool as_srgb,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 3, (const byte*)img.data(), as_srgb,
      linear, mipmap);
}
void set_texture(ogl_texture* texture, const image<vec3f>& img, bool as_float,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 3, (const float*)img.data(), as_float,
      linear, mipmap);
}

void set_texture(ogl_texture* texture, const image<byte>& img, bool as_srgb,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 1, (const byte*)img.data(), as_srgb,
      linear, mipmap);
}
void set_texture(ogl_texture* texture, const image<float>& img, bool as_float,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 1, (const float*)img.data(), as_float,
      linear, mipmap);
}

void set_cubemap(ogl_cubemap* cubemap, int size, int nchannels,
    const array<byte*, 6>& images, bool as_srgb, bool linear, bool mipmap) {
  static auto sformat = unordered_map<int, uint>{
      {1, GL_SRGB},
      {2, GL_SRGB},
      {3, GL_SRGB},
      {4, GL_SRGB_ALPHA},
  };
  static auto iformat = unordered_map<int, uint>{
      {1, GL_RGB},
      {2, GL_RGB},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  static auto cformat = unordered_map<int, uint>{
      {1, GL_RED},
      {2, GL_RG},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  assert_ogl_error();
  if (size == 0) {
    clear_cubemap(cubemap);
    return;
  }

  if (!cubemap->cubemap_id) glGenTextures(1, &cubemap->cubemap_id);
  if (cubemap->size != size || cubemap->nchannels != nchannels ||
      cubemap->is_srgb != as_srgb || cubemap->is_float == true ||
      cubemap->linear != linear || cubemap->mipmap != mipmap) {
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->cubemap_id);

    for (auto i = 0; i < 6; i++) {
      if (!images[i]) {
        throw std::runtime_error{"cannot initialize cubemap from empty image"};
        return;
      }
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0,
          as_srgb ? sformat.at(nchannels) : iformat.at(nchannels), size, size,
          0, cformat.at(nchannels), GL_UNSIGNED_BYTE, images[i]);
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER,
        linear ? GL_LINEAR : GL_NEAREST);

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    if (mipmap) {
      glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
    }
  } else {
    throw std::runtime_error{"cannot modify initialized cubemap"};
    // glBindTexture(GL_TEXTURE_2D, cubemap->cubemap_id);
    // glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
    //     cformat.at(nchannels), GL_UNSIGNED_BYTE, img);
    // if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
  cubemap->size      = size;
  cubemap->nchannels = nchannels;
  cubemap->is_srgb   = as_srgb;
  cubemap->is_float  = false;
  cubemap->linear    = linear;
  cubemap->mipmap    = mipmap;
  assert_ogl_error();
}

void set_cubemap(ogl_cubemap* cubemap, int size, int nchannels,
    const array<float*, 6>& images, bool as_float, bool linear, bool mipmap) {
  static auto fformat = unordered_map<int, uint>{
      {1, GL_RGB16F},
      {2, GL_RGB16F},
      {3, GL_RGB16F},
      {4, GL_RGBA32F},
  };
  static auto iformat = unordered_map<int, uint>{
      {1, GL_RGB},
      {2, GL_RGB},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  static auto cformat = unordered_map<int, uint>{
      {1, GL_RED},
      {2, GL_RG},
      {3, GL_RGB},
      {4, GL_RGBA},
  };
  assert_ogl_error();
  if (size == 0) {
    clear_cubemap(cubemap);
    return;
  }

  if (!cubemap->cubemap_id) glGenTextures(1, &cubemap->cubemap_id);
  if (cubemap->size != size || cubemap->nchannels != nchannels ||
      cubemap->is_float != as_float || cubemap->is_srgb == true ||
      cubemap->linear != linear || cubemap->mipmap != mipmap) {
    glGenTextures(1, &cubemap->cubemap_id);

    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->cubemap_id);

    for (auto i = 0; i < 6; i++) {
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0,
          as_float ? fformat.at(nchannels) : iformat.at(nchannels), size, size,
          0, iformat.at(nchannels), GL_FLOAT, images[i]);
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER,
        linear ? GL_LINEAR : GL_NEAREST);

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    if (mipmap) {
      glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
    }

  } else {
    // TODO(giacomo): handle this case.
    throw std::runtime_error{"cannot modify initialized cubemap"};

    //    glBindTexture(GL_TEXTURE_2D, cubemap->cubemap_id);
    //    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size, size,
    //        iformat.at(nchannels), GL_FLOAT, img);
    //    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
  cubemap->size      = size;
  cubemap->nchannels = nchannels;
  cubemap->is_srgb   = false;
  cubemap->is_float  = as_float;
  cubemap->linear    = linear;
  cubemap->mipmap    = mipmap;
  assert_ogl_error();
}

// check if cubemap is initialized
bool is_initialized(ogl_cubemap* cubemap) { return cubemap->cubemap_id != 0; }

// clear cubemap
void clear_cubemap(ogl_cubemap* cubemap) {
  if (cubemap->cubemap_id) glDeleteTextures(1, &cubemap->cubemap_id);
  *cubemap = ogl_cubemap{};
}

void set_cubemap(ogl_cubemap* cubemap, const array<image<vec4b>, 6>& img,
    int nchannels, bool as_srgb, bool linear, bool mipmap) {
  auto data = array<byte*, 6>{(byte*)img[0].data(), (byte*)img[1].data(),
      (byte*)img[2].data(), (byte*)img[3].data(), (byte*)img[4].data(),
      (byte*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_srgb, linear, mipmap);
}
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec4f>, 6>& img,
    int nchannels, bool as_float, bool linear, bool mipmap) {
  auto data = array<float*, 6>{(float*)img[0].data(), (float*)img[1].data(),
      (float*)img[2].data(), (float*)img[3].data(), (float*)img[4].data(),
      (float*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_float, linear, mipmap);
}
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec3b>, 6>& img,
    int nchannels, bool as_srgb, bool linear, bool mipmap) {
  auto data = array<byte*, 6>{(byte*)img[0].data(), (byte*)img[1].data(),
      (byte*)img[2].data(), (byte*)img[3].data(), (byte*)img[4].data(),
      (byte*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_srgb, linear, mipmap);
}
void set_cubemap(ogl_cubemap* cubemap, const array<image<vec3f>, 6>& img,
    int nchannels, bool as_float, bool linear, bool mipmap) {
  auto data = array<float*, 6>{(float*)img[0].data(), (float*)img[1].data(),
      (float*)img[2].data(), (float*)img[3].data(), (float*)img[4].data(),
      (float*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_float, linear, mipmap);
}
void set_cubemap(ogl_cubemap* cubemap, const array<image<byte>, 6>& img,
    int nchannels, bool as_srgb, bool linear, bool mipmap) {
  auto data = array<byte*, 6>{(byte*)img[0].data(), (byte*)img[1].data(),
      (byte*)img[2].data(), (byte*)img[3].data(), (byte*)img[4].data(),
      (byte*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_srgb, linear, mipmap);
}
void set_cubemap(ogl_cubemap* cubemap, const array<image<float>, 6>& img,
    int nchannels, bool as_float, bool linear, bool mipmap) {
  auto data = array<float*, 6>{(float*)img[0].data(), (float*)img[1].data(),
      (float*)img[2].data(), (float*)img[3].data(), (float*)img[4].data(),
      (float*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, nchannels, data, as_float, linear, mipmap);
}

// check if buffer is initialized
bool is_initialized(ogl_arraybuffer* buffer) { return buffer->buffer_id != 0; }

// set buffer
void set_arraybuffer(ogl_arraybuffer* buffer, size_t size, int esize,
    const float* data, bool dynamic) {
  assert_ogl_error();
  if (size == 0 || data == nullptr) {
    clear_arraybuffer(buffer);
    return;
  }
  if (!buffer->buffer_id) glGenBuffers(1, &buffer->buffer_id);
  auto target = GL_ARRAY_BUFFER;
  glBindBuffer(target, buffer->buffer_id);
  if (buffer->size != size || buffer->dynamic != dynamic) {
    glBufferData(target, size * sizeof(float), data,
        dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  } else {
    glBufferSubData(target, 0, size * sizeof(float), data);
  }
  buffer->size    = size;
  buffer->esize   = esize;
  buffer->dynamic = dynamic;
  assert_ogl_error();
}

// clear buffer
void clear_arraybuffer(ogl_arraybuffer* buffer) {
  assert_ogl_error();
  if (buffer->buffer_id) glDeleteBuffers(1, &buffer->buffer_id);
  assert_ogl_error();
  buffer->buffer_id = 0;
  buffer->size      = 0;
  buffer->esize     = 0;
  buffer->dynamic   = false;
}

// set buffer
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<float>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 1, 1, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec2f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 2, 2, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec3f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 3, 3, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec4f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 4, 4, (float*)data.data(), dynamic);
}

// set buffer
void set_elementbuffer(ogl_elementbuffer* buffer, size_t size,
    ogl_element_type element, const int* data, bool dynamic) {
  assert_ogl_error();
  if (size == 0 || data == nullptr) {
    clear_elementbuffer(buffer);
    return;
  }
  if (!buffer->buffer_id) glGenBuffers(1, &buffer->buffer_id);
  auto target = GL_ELEMENT_ARRAY_BUFFER;
  glBindBuffer(target, buffer->buffer_id);
  if (buffer->size != size || buffer->dynamic != dynamic) {
    glBufferData(target, size * sizeof(int), data,
        dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  } else {
    glBufferSubData(target, 0, size * sizeof(int), data);
  }
  buffer->size    = size;
  buffer->element = element;
  buffer->dynamic = dynamic;
  assert_ogl_error();
}

// check if buffer is initialized
bool is_initialized(ogl_elementbuffer* buffer) {
  return buffer->buffer_id != 0;
}

// clear buffer
void clear_elementbuffer(ogl_elementbuffer* buffer) {
  assert_ogl_error();
  if (buffer->buffer_id) glDeleteBuffers(1, &buffer->buffer_id);
  assert_ogl_error();
  buffer->buffer_id = 0;
  buffer->size      = 0;
  buffer->element   = ogl_element_type::points;
  buffer->dynamic   = false;
}

// set buffer
void set_elementbuffer(
    ogl_elementbuffer* buffer, const vector<int>& points, bool dynamic) {
  set_elementbuffer(buffer, points.size() * 1, ogl_element_type::points,
      (int*)points.data(), dynamic);
}
void set_elementbuffer(
    ogl_elementbuffer* buffer, const vector<vec2i>& lines, bool dynamic) {
  set_elementbuffer(buffer, lines.size() * 2, ogl_element_type::lines,
      (int*)lines.data(), dynamic);
}
void set_elementbuffer(
    ogl_elementbuffer* buffer, const vector<vec3i>& triangles, bool dynamic) {
  set_elementbuffer(buffer, triangles.size() * 3, ogl_element_type::triangles,
      (int*)triangles.data(), dynamic);
}

// initialize program
bool init_program(ogl_program* program, const string& vertex,
    const string& fragment, string& error, string& errorlog) {
  // error
  auto program_error = [&error, &errorlog, program](
                           const char* message, const char* log) {
    clear_program(program);
    error    = message;
    errorlog = log;
    return false;
  };

  // clear
  if (program->program_id) clear_program(program);

  // setup code
  program->vertex_code   = vertex;
  program->fragment_code = fragment;

  // create arrays
  assert_ogl_error();
  glGenVertexArrays(1, &program->array_id);
  glBindVertexArray(program->array_id);
  assert_ogl_error();

  const char* ccvertex   = vertex.data();
  const char* ccfragment = fragment.data();
  int         errflags;
  char        errbuf[10000];

  // create vertex
  assert_ogl_error();
  program->vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(program->vertex_id, 1, &ccvertex, NULL);
  glCompileShader(program->vertex_id);
  glGetShaderiv(program->vertex_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program->vertex_id, 10000, 0, errbuf);
    return program_error("vertex shader not compiled", errbuf);
  }
  assert_ogl_error();

  // create fragment
  assert_ogl_error();
  program->fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(program->fragment_id, 1, &ccfragment, NULL);
  glCompileShader(program->fragment_id);
  glGetShaderiv(program->fragment_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program->fragment_id, 10000, 0, errbuf);
    return program_error("fragment shader not compiled", errbuf);
  }
  assert_ogl_error();

  // create program
  assert_ogl_error();
  program->program_id = glCreateProgram();
  glAttachShader(program->program_id, program->vertex_id);
  glAttachShader(program->program_id, program->fragment_id);
  glLinkProgram(program->program_id);
  glGetProgramiv(program->program_id, GL_LINK_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program->program_id, 10000, 0, errbuf);
    return program_error("program not linked", errbuf);
  }
  // TODO(giacomo): Apparently validation must be done just before drawing.
  //    https://community.khronos.org/t/samplers-of-different-types-use-the-same-textur/66329
  // If done here, validation fails when using cubemaps and textures in the
  // same shader. We should create a function validate_program() anc call it
  // separately.
  //
  // glValidateProgram(program->program_id);
  // glGetProgramiv(program->program_id, GL_VALIDATE_STATUS, &errflags);
  // if (!errflags) {
  //   glGetProgramInfoLog(program->program_id, 10000, 0, errbuf);
  //   return program_error("program not validated", errbuf);
  // }
  assert_ogl_error();

  // done
  return true;
}

// clear program
void clear_program(ogl_program* program) {
  if (program->program_id) glDeleteProgram(program->program_id);
  if (program->vertex_id) glDeleteShader(program->vertex_id);
  if (program->fragment_id) glDeleteProgram(program->fragment_id);
  if (program->array_id) glDeleteVertexArrays(1, &program->array_id);
  program->program_id  = 0;
  program->vertex_id   = 0;
  program->fragment_id = 0;
  program->array_id    = 0;
}

bool is_initialized(const ogl_program* program) {
  return program->program_id != 0;
}

// bind program
void bind_program(ogl_program* program) {
  assert_ogl_error();
  glUseProgram(program->program_id);
  assert_ogl_error();
}
// unbind program
void unbind_program(ogl_program* program) { glUseProgram(0); }
// unbind program
void unbind_program() { glUseProgram(0); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

void init_glbuffer(
    uint& buffer_id, bool element, int size, int count, const float* array) {
  assert_ogl_error();
  glGenBuffers(1, &buffer_id);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER,
      count * size * sizeof(float), array, GL_STATIC_DRAW);
  assert_ogl_error();
}

void init_glbuffer(
    uint& buffer_id, bool element, int size, int count, const int* array) {
  assert_ogl_error();
  glGenBuffers(1, &buffer_id);
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER,
      count * size * sizeof(int), array, GL_STATIC_DRAW);
  assert_ogl_error();
}

void update_glbuffer(
    uint& buffer_id, bool element, int size, int count, const int* array) {
  assert_ogl_error();
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferSubData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, 0,
      size * count * sizeof(int), array);
  assert_ogl_error();
}

void update_glbuffer(
    uint& buffer_id, bool element, int size, int count, const float* array) {
  assert_ogl_error();
  glBindBuffer(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, buffer_id);
  glBufferSubData(element ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER, 0,
      size * count * sizeof(float), array);
  assert_ogl_error();
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

ogl_image::~ogl_image() {
  if (program) delete program;
  if (texcoords) delete texcoords;
  if (triangles) delete triangles;
}

bool is_initialized(const ogl_image* image) {
  return is_initialized(image->program);
}

// init image program
bool init_image(ogl_image* image) {
  if (is_initialized(image)) return true;

  auto texcoords = vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  auto triangles = vector<vec3i>{{0, 1, 2}, {0, 2, 3}};

  auto error = ""s, errorlog = ""s;
  if (!init_program(
          image->program, glimage_vertex, glimage_fragment, error, errorlog))
    return false;
  set_arraybuffer(
      image->texcoords, texcoords.size() * 2, 2, (float*)texcoords.data());
  set_elementbuffer(image->triangles, triangles.size() * 3,
      ogl_element_type::triangles, (int*)triangles.data());
  return true;
}

// clear an opengl image
void clear_image(ogl_image* image) {
  clear_program(image->program);
  clear_texture(image->texture);
  clear_arraybuffer(image->texcoords);
  clear_elementbuffer(image->triangles);
}

// update image data
void set_image(
    ogl_image* oimg, const image<vec4f>& img, bool linear, bool mipmap) {
  set_texture(oimg->texture, img, false, linear, mipmap);
}
void set_image(
    ogl_image* oimg, const image<vec4b>& img, bool linear, bool mipmap) {
  set_texture(oimg->texture, img, false, linear, mipmap);
}

// draw image
void draw_image(ogl_image* image, const ogl_image_params& params) {
  assert_ogl_error();
  glViewport(params.framebuffer.x, params.framebuffer.y, params.framebuffer.z,
      params.framebuffer.w);
  glClearColor(params.background.x, params.background.y, params.background.z,
      params.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  bind_program(image->program);
  set_uniform(image->program, "txt", image->texture, 0);
  set_uniform(image->program, "window_size",
      vec2f{(float)params.window.x, (float)params.window.y});
  set_uniform(image->program, "image_size",
      vec2f{(float)image->texture->size.x, (float)image->texture->size.y});
  set_uniform(image->program, "image_center", params.center);
  set_uniform(image->program, "image_scale", params.scale);
  set_attribute(image->program, "texcoord", image->texcoords);
  draw_elements(image->triangles);
  unbind_program(image->program);
  assert_ogl_error();
}

// set uniforms
void set_uniform(ogl_program* program, int location, int value) {
  assert_ogl_error();
  glUniform1i(location, value);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec2i& value) {
  assert_ogl_error();
  glUniform2i(location, value.x, value.y);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec3i& value) {
  assert_ogl_error();
  glUniform3i(location, value.x, value.y, value.z);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec4i& value) {
  assert_ogl_error();
  glUniform4i(location, value.x, value.y, value.z, value.w);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, float value) {
  assert_ogl_error();
  glUniform1f(location, value);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec2f& value) {
  assert_ogl_error();
  glUniform2f(location, value.x, value.y);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec3f& value) {
  assert_ogl_error();
  glUniform3f(location, value.x, value.y, value.z);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const vec4f& value) {
  assert_ogl_error();
  glUniform4f(location, value.x, value.y, value.z, value.w);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const mat2f& value) {
  assert_ogl_error();
  glUniformMatrix2fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const mat3f& value) {
  assert_ogl_error();
  glUniformMatrix3fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const mat4f& value) {
  assert_ogl_error();
  glUniformMatrix4fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const frame2f& value) {
  assert_ogl_error();
  glUniformMatrix3x2fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(ogl_program* program, int location, const frame3f& value) {
  assert_ogl_error();
  glUniformMatrix4x3fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

// get uniform location
int get_uniform_location(ogl_program* program, const char* name) {
  return glGetUniformLocation(program->program_id, name);
}

// set uniform texture
void set_uniform(
    ogl_program* program, int location, const ogl_texture* texture, int unit) {
  assert_ogl_error();
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture->texture_id);
  glUniform1i(location, unit);
  assert_ogl_error();
}
void set_uniform(ogl_program* program, const char* name,
    const ogl_texture* texture, int unit) {
  return set_uniform(
      program, get_uniform_location(program, name), texture, unit);
}
void set_uniform(ogl_program* program, int location, int location_on,
    const ogl_texture* texture, int unit) {
  assert_ogl_error();
  if (texture && texture->texture_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture->texture_id);
    glUniform1i(location, unit);
    glUniform1i(location_on, 1);
  } else {
    glUniform1i(location_on, 0);
  }
  assert_ogl_error();
}
void set_uniform(ogl_program* program, const char* name, const char* name_on,
    const ogl_texture* texture, int unit) {
  return set_uniform(program, get_uniform_location(program, name),
      get_uniform_location(program, name_on), texture, unit);
}

// set uniform cubemap
void set_uniform(
    ogl_program* program, int location, const ogl_cubemap* cubemap, int unit) {
  assert_ogl_error();
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->cubemap_id);
  glUniform1i(location, unit);
  assert_ogl_error();
}
void set_uniform(ogl_program* program, const char* name,
    const ogl_cubemap* cubemap, int unit) {
  return set_uniform(
      program, get_uniform_location(program, name), cubemap, unit);
}
void set_uniform(ogl_program* program, int location, int location_on,
    const ogl_cubemap* cubemap, int unit) {
  assert_ogl_error();
  if (cubemap && cubemap->cubemap_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->cubemap_id);
    glUniform1i(location, unit);
    glUniform1i(location_on, 1);
  } else {
    glUniform1i(location_on, 0);
  }
  assert_ogl_error();
}
void set_uniform(ogl_program* program, const char* name, const char* name_on,
    const ogl_cubemap* cubemap, int unit) {
  return set_uniform(program, get_uniform_location(program, name),
      get_uniform_location(program, name_on), cubemap, unit);
}

// get attribute location
int get_attribute_location(ogl_program* program, const char* name) {
  return glGetAttribLocation(program->program_id, name);
}

// set vertex attributes
void set_attribute(
    ogl_program* program, int location, ogl_arraybuffer* buffer) {
  assert_ogl_error();
  glBindBuffer(GL_ARRAY_BUFFER, buffer->buffer_id);
  glEnableVertexAttribArray(location);
  glVertexAttribPointer(location, buffer->esize, GL_FLOAT, false, 0, nullptr);
  assert_ogl_error();
}
void set_attribute(
    ogl_program* program, const char* name, ogl_arraybuffer* buffer) {
  return set_attribute(program, get_attribute_location(program, name), buffer);
}

// set vertex attributes
void set_attribute(ogl_program* program, int location, float value) {
  glDisableVertexAttribArray(location);
  glVertexAttrib1f(location, value);
}
void set_attribute(ogl_program* program, int location, const vec2f& value) {
  glDisableVertexAttribArray(location);
  glVertexAttrib2f(location, value.x, value.y);
}
void set_attribute(ogl_program* program, int location, const vec3f& value) {
  glDisableVertexAttribArray(location);
  glVertexAttrib3f(location, value.x, value.y, value.z);
}
void set_attribute(ogl_program* program, int location, const vec4f& value) {
  glDisableVertexAttribArray(location);
  glVertexAttrib4f(location, value.x, value.y, value.z, value.w);
}

// draw elements
void draw_elements(ogl_elementbuffer* buffer) {
  static auto elements = unordered_map<ogl_element_type, uint>{
      {ogl_element_type::points, GL_POINTS},
      {ogl_element_type::lines, GL_LINES},
      {ogl_element_type::triangles, GL_TRIANGLES},
  };
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer->buffer_id);
  glDrawElements(elements.at(buffer->element), (GLsizei)buffer->size,
      GL_UNSIGNED_INT, nullptr);
}

void set_framebuffer(ogl_framebuffer* framebuffer, const vec2i& size) {
  assert(!framebuffer->framebuffer_id);
  assert(!framebuffer->renderbuffer_id);

  glGenFramebuffers(1, &framebuffer->framebuffer_id);
  glGenRenderbuffers(1, &framebuffer->renderbuffer_id);
  framebuffer->size = size;

  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer->framebuffer_id);
  glBindRenderbuffer(GL_RENDERBUFFER, framebuffer->renderbuffer_id);

  // create render buffer for depth and stencil
  // TODO(giacomo): Why do we need to put STENCIL8 to make things work on Mac??
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8,
      framebuffer->size.x, framebuffer->size.y);

  // bind together frame buffer and render buffer
  // TODO(giacomo): We put STENCIL here for the same reason...
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT,
      GL_RENDERBUFFER, framebuffer->renderbuffer_id);
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);

  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  assert_ogl_error();
}

inline void set_framebuffer_texture(const ogl_framebuffer* framebuffer,
    uint texture_id, uint target, uint mipmap_level) {
  // TODO(giacomo): We change the state of the framebuffer, but we don't store
  // this information anywhere, unlike the rest of the library.
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer->framebuffer_id);
  glFramebufferTexture2D(
      GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, target, texture_id, mipmap_level);
  // glBindFramebuffer(GL_FRAMEBUFFER, ogl_framebuffer::bound_framebuffer_id);
  assert_ogl_error();
}

bool is_framebuffer_bound(const ogl_framebuffer* framebuffer) {
  return framebuffer->framebuffer_id == ogl_framebuffer::bound_framebuffer_id;
}

void set_framebuffer_texture(const ogl_framebuffer* framebuffer,
    const ogl_texture* texture, uint mipmap_level) {
  set_framebuffer_texture(
      framebuffer, texture->texture_id, GL_TEXTURE_2D, mipmap_level);
}

void set_framebuffer_texture(const ogl_framebuffer* framebuffer,
    const ogl_cubemap* cubemap, uint face, uint mipmap_level) {
  set_framebuffer_texture(framebuffer, cubemap->cubemap_id,
      GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, mipmap_level);
}

void bind_framebuffer(const ogl_framebuffer* framebuffer) {
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer->framebuffer_id);
  ogl_framebuffer::bound_framebuffer_id = framebuffer->framebuffer_id;
  assert_ogl_error();
}

void unbind_framebuffer() {
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  ogl_framebuffer::bound_framebuffer_id = 0;
  assert_ogl_error();
}

void clear_framebuffer(ogl_framebuffer* framebuffer) {
  glDeleteFramebuffers(1, &framebuffer->framebuffer_id);
  glDeleteRenderbuffers(1, &framebuffer->renderbuffer_id);
  *framebuffer = {};
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

static const char* glscene_fragment =
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
  if(mtype == 0) {
    ke = emission;
    kd = vec3(0,0,0);
    ks = vec3(0,0,0);
    op = 1;
    return false;
  }

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

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

// forward declaration
void clear_shape(ogl_shape* shape);

ogl_shape::~ogl_shape() {
  clear_shape(this);
  if (positions) delete positions;
  if (normals) delete normals;
  if (texcoords) delete texcoords;
  if (colors) delete colors;
  if (tangents) delete tangents;
  if (points) delete points;
  if (lines) delete lines;
  if (triangles) delete triangles;
  if (quads) delete quads;
  if (edges) delete edges;
}

ogl_scene::~ogl_scene() {
  clear_scene(this);
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto texture : textures) delete texture;
  for (auto light : lights) delete light;
}

// Initialize an OpenGL scene
void init_scene(ogl_scene* scene) {
  if (is_initialized(scene->program)) return;
  auto error = ""s, errorlog = ""s;
  init_program(
      scene->program, glscene_vertex, glscene_fragment, error, errorlog);
}
bool is_initialized(ogl_scene* scene) { return is_initialized(scene->program); }

// Clear an OpenGL shape
void clear_shape(ogl_shape* shape) {
  clear_arraybuffer(shape->positions);
  clear_arraybuffer(shape->normals);
  clear_arraybuffer(shape->texcoords);
  clear_arraybuffer(shape->colors);
  clear_arraybuffer(shape->tangents);
  clear_elementbuffer(shape->points);
  clear_elementbuffer(shape->lines);
  clear_elementbuffer(shape->triangles);
  clear_elementbuffer(shape->quads);
  clear_elementbuffer(shape->edges);
}

// Clear an OpenGL scene
void clear_scene(ogl_scene* scene) {
  for (auto texture : scene->textures) clear_texture(texture);
  for (auto shape : scene->shapes) clear_shape(shape);
  clear_program(scene->program);
  clear_program(scene->environment_program);
  clear_cubemap(scene->environment_cubemap);
}

// add camera
ogl_camera* add_camera(ogl_scene* scene) {
  return scene->cameras.emplace_back(new ogl_camera{});
}
void set_frame(ogl_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(ogl_camera* camera, float lens, float aspect, float film) {
  camera->lens   = lens;
  camera->aspect = aspect;
  camera->film   = film;
}
void set_nearfar(ogl_camera* camera, float near, float far) {
  camera->near = near;
  camera->far  = far;
}

// add texture
ogl_texture* add_texture(ogl_scene* scene) {
  return scene->textures.emplace_back(new ogl_texture{});
}

// add shape
ogl_shape* add_shape(ogl_scene* scene) {
  return scene->shapes.emplace_back(new ogl_shape{});
}

void set_points(ogl_shape* shape, const vector<int>& points) {
  set_elementbuffer(shape->points, points);
}
void set_lines(ogl_shape* shape, const vector<vec2i>& lines) {
  set_elementbuffer(shape->lines, lines);
}
void set_triangles(ogl_shape* shape, const vector<vec3i>& triangles) {
  set_elementbuffer(shape->triangles, triangles);
}
void set_quads(ogl_shape* shape, const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  set_elementbuffer(shape->quads, triangles);
}
void set_edges(ogl_shape* shape, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
  auto edgemap = unordered_set<vec2i>{};
  for (auto t : triangles) {
    edgemap.insert({min(t.x, t.y), max(t.x, t.y)});
    edgemap.insert({min(t.y, t.z), max(t.y, t.z)});
    edgemap.insert({min(t.z, t.x), max(t.z, t.x)});
  }
  for (auto t : quads) {
    edgemap.insert({min(t.x, t.y), max(t.x, t.y)});
    edgemap.insert({min(t.y, t.z), max(t.y, t.z)});
    edgemap.insert({min(t.z, t.w), max(t.z, t.w)});
    edgemap.insert({min(t.w, t.x), max(t.w, t.x)});
  }
  auto edges = vector<vec2i>(edgemap.begin(), edgemap.end());
  set_elementbuffer(shape->edges, edges);
}
void set_positions(ogl_shape* shape, const vector<vec3f>& positions) {
  set_arraybuffer(shape->positions, positions);
}
void set_normals(ogl_shape* shape, const vector<vec3f>& normals) {
  set_arraybuffer(shape->normals, normals);
}
void set_texcoords(ogl_shape* shape, const vector<vec2f>& texcoords) {
  set_arraybuffer(shape->texcoords, texcoords);
}
void set_colors(ogl_shape* shape, const vector<vec3f>& colors) {
  set_arraybuffer(shape->colors, colors);
}
void set_colors(ogl_shape* shape, const vector<vec4f>& colors) {
  set_arraybuffer(shape->colors, colors);
}
void set_tangents(ogl_shape* shape, const vector<vec4f>& tangents) {
  set_arraybuffer(shape->tangents, tangents);
}

// add instance
ogl_instance* add_instance(ogl_scene* scene) {
  return scene->instances.emplace_back(new ogl_instance{});
}
void set_frame(ogl_instance* instance, const frame3f& frame) {
  instance->frame = frame;
}
void set_shape(ogl_instance* instance, ogl_shape* shape) {
  instance->shape = shape;
}
void set_material(ogl_instance* instance, ogl_material* material) {
  instance->material = material;
}
void set_hidden(ogl_instance* instance, bool hidden) {
  instance->hidden = hidden;
}
void set_highlighted(ogl_instance* instance, bool highlighted) {
  instance->highlighted = highlighted;
}

// add material
ogl_material* add_material(ogl_scene* scene) {
  return scene->materials.emplace_back(new ogl_material{});
}
void set_emission(
    ogl_material* material, const vec3f& emission, ogl_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    ogl_material* material, const vec3f& color, ogl_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    ogl_material* material, float specular, ogl_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_roughness(
    ogl_material* material, float roughness, ogl_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    ogl_material* material, float opacity, ogl_texture* opacity_tex) {
  material->opacity = opacity;
}
void set_metallic(
    ogl_material* material, float metallic, ogl_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_normalmap(ogl_material* material, ogl_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// shortcuts
ogl_camera* add_camera(ogl_scene* scene, const frame3f& frame, float lens,
    float aspect, float film, float near, float far) {
  auto camera = add_camera(scene);
  set_frame(camera, frame);
  set_lens(camera, lens, aspect, film);
  set_nearfar(camera, near, far);
  return camera;
}
ogl_material* add_material(ogl_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    ogl_texture* emission_tex, ogl_texture* color_tex,
    ogl_texture* specular_tex, ogl_texture* metallic_tex,
    ogl_texture* roughness_tex, ogl_texture* normalmap_tex) {
  auto material = add_material(scene);
  set_emission(material, emission, emission_tex);
  set_color(material, color, color_tex);
  set_specular(material, specular, specular_tex);
  set_metallic(material, metallic, metallic_tex);
  set_roughness(material, roughness, roughness_tex);
  set_normalmap(material, normalmap_tex);
  return material;
}
ogl_shape* add_shape(ogl_scene* scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec4f>& colors, bool edges) {
  auto shape = add_shape(scene);
  set_points(shape, points);
  set_lines(shape, lines);
  set_triangles(shape, triangles);
  set_quads(shape, quads);
  set_positions(shape, positions);
  set_normals(shape, normals);
  set_texcoords(shape, texcoords);
  set_colors(shape, colors);
  if (edges && (!triangles.empty() || !quads.empty())) {
    set_edges(shape, triangles, quads);
  }
  return shape;
}
ogl_instance* add_instance(ogl_scene* scene, const frame3f& frame,
    ogl_shape* shape, ogl_material* material, bool hidden, bool highlighted) {
  auto instance = add_instance(scene);
  set_frame(instance, frame);
  set_shape(instance, shape);
  set_material(instance, material);
  set_hidden(instance, hidden);
  set_highlighted(instance, highlighted);
  return instance;
}

// add light
ogl_light* add_light(ogl_scene* scene) {
  return scene->lights.emplace_back(new ogl_light{});
}
void set_light(ogl_light* light, const vec3f& position, const vec3f& emission,
    ogl_light_type type, bool camera) {
  light->position = position;
  light->emission = emission;
  light->type     = type;
  light->camera   = camera;
}
void clear_lights(ogl_scene* scene) {
  for (auto light : scene->lights) delete light;
  scene->lights.clear();
}
bool has_max_lights(ogl_scene* scene) { return scene->lights.size() >= 16; }
void add_default_lights(ogl_scene* scene) {
  clear_lights(scene);
  set_light(add_light(scene), normalize(vec3f{1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{-1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{-1, -1, 1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{0.1, 0.5, -1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true);
}

// Draw a shape
void draw_object(
    ogl_scene* scene, ogl_instance* instance, const ogl_scene_params& params) {
  static auto empty_instances = vector<frame3f>{identity3x4f};

  if (instance->hidden) return;

  assert_ogl_error();
  auto shape_xform     = frame_to_mat(instance->frame);
  auto shape_inv_xform = transpose(
      frame_to_mat(inverse(instance->frame, params.non_rigid_frames)));
  set_uniform(scene->program, "frame", shape_xform);
  set_uniform(scene->program, "frameit", shape_inv_xform);
  set_uniform(scene->program, "offset", 0.0f);
  if (instance->highlighted) {
    set_uniform(scene->program, "highlight", vec4f{1, 1, 0, 1});
  } else {
    set_uniform(scene->program, "highlight", vec4f{0, 0, 0, 0});
  }
  assert_ogl_error();

  auto material = instance->material;
  auto mtype    = 2;
  set_uniform(scene->program, "mtype", mtype);
  set_uniform(scene->program, "emission", material->emission);
  set_uniform(scene->program, "diffuse", material->color);
  set_uniform(scene->program, "specular",
      vec3f{material->metallic, material->metallic, material->metallic});
  set_uniform(scene->program, "roughness", material->roughness);
  set_uniform(scene->program, "opacity", material->opacity);
  set_uniform(scene->program, "double_sided", (int)params.double_sided);
  set_uniform(scene->program, "emission_tex", "emission_tex_on",
      material->emission_tex, 0);
  set_uniform(
      scene->program, "diffuse_tex", "diffuse_tex_on", material->color_tex, 1);
  set_uniform(scene->program, "specular_tex", "specular_tex_on",
      material->metallic_tex, 2);
  set_uniform(scene->program, "roughness_tex", "roughness_tex_on",
      material->roughness_tex, 3);
  set_uniform(scene->program, "opacity_tex", "opacity_tex_on",
      material->opacity_tex, 4);
  set_uniform(scene->program, "mat_norm_tex", "mat_norm_tex_on",
      material->normal_tex, 5);
  assert_ogl_error();

  auto shape = instance->shape;
  set_uniform(scene->program, "faceted", !is_initialized(shape->normals));
  set_attribute(scene->program, "positions", shape->positions, vec3f{0, 0, 0});
  set_attribute(scene->program, "normals", shape->normals, vec3f{0, 0, 1});
  set_attribute(scene->program, "texcoords", shape->texcoords, vec2f{0, 0});
  set_attribute(scene->program, "colors", shape->colors, vec4f{1, 1, 1, 1});
  set_attribute(scene->program, "tangents", shape->tangents, vec4f{0, 0, 1, 1});
  assert_ogl_error();

  if (is_initialized(shape->points)) {
    glPointSize(shape->points_size);
    set_uniform(scene->program, "etype", 1);
    draw_elements(shape->points);
  }
  if (is_initialized(shape->lines)) {
    set_uniform(scene->program, "etype", 2);
    draw_elements(shape->lines);
  }
  if (is_initialized(shape->triangles)) {
    set_uniform(scene->program, "etype", 3);
    draw_elements(shape->triangles);
  }
  if (is_initialized(shape->quads)) {
    set_uniform(scene->program, "etype", 3);
    draw_elements(shape->quads);
  }
  assert_ogl_error();

  if (is_initialized(shape->edges) && params.edges && !params.wireframe) {
    set_uniform(scene->program, "mtype", mtype);
    set_uniform(scene->program, "emission", vec3f{0, 0, 0});
    set_uniform(scene->program, "diffuse", vec3f{0, 0, 0});
    set_uniform(scene->program, "specular", vec3f{0, 0, 0});
    set_uniform(scene->program, "roughness", 0);
    set_uniform(scene->program, "etype", 3);
    draw_elements(shape->edges);
    assert_ogl_error();
  }
}

// Display a scene
void draw_scene(ogl_scene* scene, ogl_camera* camera, const vec4i& viewport,
    const ogl_scene_params& params) {
  static auto camera_light0 = ogl_light{normalize(vec3f{1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true};
  static auto camera_light1 = ogl_light{normalize(vec3f{-1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true};
  static auto camera_light2 = ogl_light{normalize(vec3f{-1, -1, 1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true};
  static auto camera_light3 = ogl_light{normalize(vec3f{0.1, 0.5, -1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true};
  static auto camera_lights = vector<ogl_light*>{
      &camera_light0, &camera_light1, &camera_light2, &camera_light3};
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * atan(camera->film / (camera_aspect * 2 * camera->lens)))
          : (2 * atan(camera->film / (2 * camera->lens)));
  auto camera_view = frame_to_mat(inverse(camera->frame));
  auto camera_proj = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  assert_ogl_error();
  clear_ogl_framebuffer(params.background);
  set_ogl_viewport(viewport);

  assert_ogl_error();
  bind_program(scene->program);
  assert_ogl_error();
  set_uniform(scene->program, "eye", camera->frame.o);
  set_uniform(scene->program, "view", camera_view);
  set_uniform(scene->program, "projection", camera_proj);
  set_uniform(scene->program, "eyelight",
      params.shading == ogl_shading_type::eyelight ? 1 : 0);
  set_uniform(scene->program, "exposure", params.exposure);
  set_uniform(scene->program, "gamma", params.gamma);
  assert_ogl_error();

  if (params.shading == ogl_shading_type::lights ||
      params.shading == ogl_shading_type::camlights) {
    assert_ogl_error();
    auto& lights = params.shading == ogl_shading_type::lights ? scene->lights
                                                              : camera_lights;
    set_uniform(scene->program, "lamb", vec3f{0, 0, 0});
    set_uniform(scene->program, "lnum", (int)lights.size());
    auto lid = 0;
    for (auto light : lights) {
      auto is = std::to_string(lid);
      if (light->camera) {
        auto position = light->type == ogl_light_type::directional
                            ? transform_direction(
                                  camera->frame, light->position)
                            : transform_point(camera->frame, light->position);
        set_uniform(scene->program, ("lpos[" + is + "]").c_str(), position);
      } else {
        set_uniform(
            scene->program, ("lpos[" + is + "]").c_str(), light->position);
      }
      set_uniform(scene->program, ("lke[" + is + "]").c_str(), light->emission);
      set_uniform(
          scene->program, ("ltype[" + is + "]").c_str(), (int)light->type);
      lid++;
    }
    assert_ogl_error();
  }

  if (params.wireframe) set_ogl_wireframe(true);
  for (auto instance : scene->instances) {
    draw_object(scene, instance, params);
  }

  unbind_program();
  if (params.wireframe) set_ogl_wireframe(false);
}

}  // namespace yocto
