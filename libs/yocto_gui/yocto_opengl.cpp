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

#include <yocto/yocto_commonio.h>

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

GLenum _assert_ogl_error() {
  auto error_code = glGetError();
  if (error_code != GL_NO_ERROR) {
    auto error = ""s;
    switch (error_code) {
      case GL_INVALID_ENUM: error = "INVALID_ENUM"; break;
      case GL_INVALID_VALUE: error = "INVALID_VALUE"; break;
      case GL_INVALID_OPERATION: error = "INVALID_OPERATION"; break;
      // case GL_STACK_OVERFLOW: error = "STACK_OVERFLOW"; break;
      // case GL_STACK_UNDERFLOW: error = "STACK_UNDERFLOW"; break;
      case GL_OUT_OF_MEMORY: error = "OUT_OF_MEMORY"; break;
      case GL_INVALID_FRAMEBUFFER_OPERATION:
        error = "INVALID_FRAMEBUFFER_OPERATION";
        break;
    }
    printf("\n    OPENGL ERROR: %s\n\n", error.c_str());
  }
  return error_code;
}
void assert_ogl_error() { assert(_assert_ogl_error() == GL_NO_ERROR); }

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
  program->program_id  = 0;
  program->vertex_id   = 0;
  program->fragment_id = 0;
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
  assert_ogl_error();
}

void set_framebuffer(ogl_framebuffer* framebuffer, const vec2i& size) {
  if (!framebuffer->framebuffer_id) {
    glGenFramebuffers(1, &framebuffer->framebuffer_id);
  }

  if (!framebuffer->renderbuffer_id) {
    glGenRenderbuffers(1, &framebuffer->renderbuffer_id);
    // bind together frame buffer and render buffer
    // TODO(giacomo): We put STENCIL here for the same reason...
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer->framebuffer_id);
    glBindRenderbuffer(GL_RENDERBUFFER, framebuffer->renderbuffer_id);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT,
        GL_RENDERBUFFER, framebuffer->renderbuffer_id);
  }

  if (size != framebuffer->size) {
    // create render buffer for depth and stencil
    // TODO(giacomo): Why do we need to put STENCIL8 to make things work on
    // Mac??
    glBindRenderbuffer(GL_RENDERBUFFER, framebuffer->renderbuffer_id);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, size.x, size.y);
    framebuffer->size = size;
  }

  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer->framebuffer_id);
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
  glBindFramebuffer(GL_FRAMEBUFFER, ogl_framebuffer::bound_framebuffer_id);

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

void bind_shape(const ogl_shape* shape) { glBindVertexArray(shape->shape_id); }

void set_shape(ogl_shape* shape) {
  glGenVertexArrays(1, &shape->shape_id);
  assert_ogl_error();
}

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
  glDeleteVertexArrays(1, &shape->shape_id);
  shape->shape_id = 0;
}

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

void set_vertex_attribute(int location, float value) {
  glVertexAttrib1f(location, value);
}
void set_vertex_attribute(int location, const vec2f& value) {
  glVertexAttrib2f(location, value.x, value.y);
}
void set_vertex_attribute(int location, const vec3f& value) {
  glVertexAttrib3f(location, value.x, value.y, value.z);
}
void set_vertex_attribute(int location, const vec4f& value) {
  glVertexAttrib4f(location, value.x, value.y, value.z, value.w);
}

void set_vertex_attribute(int location, const ogl_arraybuffer* buffer) {
  assert_ogl_error();
  assert(buffer->buffer_id);  // == 0) return;
  glBindBuffer(GL_ARRAY_BUFFER, buffer->buffer_id);
  glEnableVertexAttribArray(location);
  glVertexAttribPointer(location, buffer->esize, GL_FLOAT, false, 0, nullptr);
  assert_ogl_error();
}

template <typename T>
void set_vertex_attribute(ogl_shape* shape, ogl_arraybuffer* attribute,
    const vector<T>& data, int location) {
  assert(!data.empty());
  set_arraybuffer(attribute, data, false);
  assert_ogl_error();
  bind_shape(shape);
  set_vertex_attribute(location, attribute);
  assert_ogl_error();
}

template <typename T>
void set_vertex_attribute(ogl_shape* shape, const T& attribute, int location) {
  assert_ogl_error();
  bind_shape(shape);
  set_vertex_attribute(location, attribute);
  assert_ogl_error();
}

void set_positions(ogl_shape* shape, const vector<vec3f>& positions) {
  if (positions.empty())
    set_vertex_attribute(shape, vec3f{0, 0, 0}, 0);
  else
    set_vertex_attribute(shape, shape->positions, positions, 0);
}
void set_normals(ogl_shape* shape, const vector<vec3f>& normals) {
  if (normals.empty())
    set_vertex_attribute(shape, vec3f{0, 0, 1}, 1);
  else
    set_vertex_attribute(shape, shape->normals, normals, 1);
}
void set_texcoords(ogl_shape* shape, const vector<vec2f>& texcoords) {
  if (texcoords.empty())
    set_vertex_attribute(shape, vec2f{0, 0}, 2);
  else
    set_vertex_attribute(shape, shape->texcoords, texcoords, 2);
}
void set_colors(ogl_shape* shape, const vector<vec4f>& colors) {
  if (colors.empty())
    set_vertex_attribute(shape, vec4f{1, 1, 1, 1}, 3);
  else
    set_vertex_attribute(shape, shape->colors, colors, 3);
}
void set_tangents(ogl_shape* shape, const vector<vec4f>& tangents) {
  if (tangents.empty())
    set_vertex_attribute(shape, vec4f{0, 0, 1, 1}, 4);
  else
    set_vertex_attribute(shape, shape->tangents, tangents, 4);
}

template <typename T>
void set_index_buffer(
    ogl_shape* shape, ogl_elementbuffer* elements, const vector<T>& data) {
  bind_shape(shape);
  set_arraybuffer(elements, data);
  assert_ogl_error();
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

void draw_shape(const ogl_shape* shape) {
  bind_shape(shape);

  if (is_initialized(shape->points)) {
    glPointSize(shape->points_size);
    draw_elements(shape->points);
  } else if (is_initialized(shape->lines)) {
    draw_elements(shape->lines);
  } else if (is_initialized(shape->triangles)) {
    draw_elements(shape->triangles);
  } else if (is_initialized(shape->quads)) {
    draw_elements(shape->quads);
  }
}

ogl_shape* cube_shape() {
  static ogl_shape* cube = nullptr;
  if (!cube) {
    // clang-format off
    static const auto cube_positions = vector<vec3f>{
      {1, -1, -1}, {1, -1,  1}, {-1, -1,  1}, {-1, -1, -1},
      {1,  1, -1}, {1,  1,  1}, {-1,  1,  1}, {-1,  1, -1},
    };
    static const auto cube_triangles = vector<vec3i>{
      {1, 3, 0}, {7, 5, 4}, {4, 1, 0}, {5, 2, 1},
      {2, 7, 3}, {0, 7, 4}, {1, 2, 3}, {7, 6, 5},
      {4, 5, 1}, {5, 6, 2}, {2, 6, 7}, {0, 3, 7}
    };
    // clang-format on
    cube = new ogl_shape{};
    set_shape(cube);
    set_positions(cube, cube_positions);
    set_triangles(cube, cube_triangles);
  }
  return cube;
}

ogl_shape* quad_shape() {
  static ogl_shape* quad = nullptr;
  if (!quad) {
    // clang-format off
    static const auto quad_positions = vector<vec3f>{
      {-1, -1, 0}, {1, -1,  0}, {1, 1,  0}, {-1, 1, 0},
    };
    static const auto quad_triangles = vector<vec3i>{
      {0, 1, 3}, {3, 2, 1}
    };
    // clang-format on
    quad = new ogl_shape{};
    set_shape(quad);
    set_positions(quad, quad_positions);
    set_triangles(quad, quad_triangles);
  }
  return quad;
}

ogl_image::~ogl_image() {
  if (program) delete program;
  if (quad) delete quad;
}

bool is_initialized(const ogl_image* image) {
  return is_initialized(image->program);
}

// init image program
bool init_image(ogl_image* image) {
  if (is_initialized(image)) return true;
  auto error = ""s, errorlog = ""s;
  if (!init_program(
          image->program, glimage_vertex, glimage_fragment, error, errorlog))
    return false;

  auto texcoords = vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  auto triangles = vector<vec3i>{{0, 1, 2}, {0, 2, 3}};
  set_shape(image->quad);
  set_vertex_attribute(image->quad, image->quad->texcoords, texcoords, 0);
  set_triangles(image->quad, triangles);

  return true;
}

// clear an opengl image
void clear_image(ogl_image* image) {
  clear_program(image->program);
  clear_texture(image->texture);
  clear_shape(image->quad);
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
  draw_shape(image->quad);
  unbind_program(image->program);
  assert_ogl_error();
}

}  // namespace yocto
