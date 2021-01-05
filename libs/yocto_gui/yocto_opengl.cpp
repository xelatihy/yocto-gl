//
// Utilities to use OpenGL 3, GLFW and ImGui.
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

#include "yocto_opengl.h"

#include <cassert>
#include <stdexcept>
#include <unordered_set>

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
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTOR HASHING
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::vec2i> {
  size_t operator()(const yocto::vec2i& v) const {
    static const auto hasher = std::hash<int>();
    auto              h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

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
void set_ogl_msaa() { glEnable(GL_MULTISAMPLE); }

void set_texture(ogl_texture& texture, const vec2i& size, int num_channels,
    const byte* img, bool as_srgb, bool linear, bool mipmap, bool wrap_repeat) {
  static auto sformat = vector<uint>{
      0, GL_SRGB, GL_SRGB, GL_SRGB, GL_SRGB_ALPHA};
  static auto iformat = vector<uint>{0, GL_RGB, GL_RGB, GL_RGB, GL_RGBA};
  static auto cformat = vector<uint>{0, GL_RED, GL_RG, GL_RGB, GL_RGBA};
  assert_ogl_error();

  if (size == zero2i) {
    clear_texture(texture);
    return;
  }

  if (!texture.texture_id) glGenTextures(1, &texture.texture_id);
  if (texture.size != size || texture.num_channels != num_channels ||
      texture.is_srgb != as_srgb || texture.is_float == true ||
      texture.linear != linear || texture.mipmap != mipmap) {
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0,
        as_srgb ? sformat.at(num_channels) : iformat.at(num_channels), size.x,
        size.y, 0, cformat.at(num_channels), GL_UNSIGNED_BYTE, img);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
    if (mipmap && img) glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
        cformat.at(num_channels), GL_UNSIGNED_BYTE, img);
    if (mipmap && img) glGenerateMipmap(GL_TEXTURE_2D);
  }
  if (wrap_repeat) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  }
  texture.size         = size;
  texture.num_channels = num_channels;
  texture.is_srgb      = as_srgb;
  texture.is_float     = false;
  texture.linear       = linear;
  texture.mipmap       = mipmap;
  assert_ogl_error();
}

void set_texture(ogl_texture& texture, const vec2i& size, int num_channels,
    const float* img, bool as_float, bool linear, bool mipmap,
    bool wrap_repeat) {
  static auto fformat = vector<uint>{
      0, GL_RGB16F, GL_RGB16F, GL_RGB16F, GL_RGBA32F};
  static auto iformat = vector<uint>{0, GL_RGB, GL_RGB, GL_RGB, GL_RGBA};
  static auto cformat = vector<uint>{0, GL_RED, GL_RG, GL_RGB, GL_RGBA};
  assert_ogl_error();

  if (size == zero2i) {
    clear_texture(texture);
    return;
  }

  if (!texture.texture_id) glGenTextures(1, &texture.texture_id);
  if (texture.size != size || texture.num_channels != num_channels ||
      texture.is_float != as_float || texture.is_srgb == true ||
      texture.linear != linear || texture.mipmap != mipmap) {
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexImage2D(GL_TEXTURE_2D, 0,
        as_float ? fformat.at(num_channels) : iformat.at(num_channels), size.x,
        size.y, 0, iformat.at(num_channels), GL_FLOAT, img);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        mipmap ? (linear ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST)
               : (linear ? GL_LINEAR : GL_NEAREST));
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
    if (mipmap && img) glGenerateMipmap(GL_TEXTURE_2D);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
        iformat.at(num_channels), GL_FLOAT, img);
    if (mipmap && img) glGenerateMipmap(GL_TEXTURE_2D);
  }
  if (wrap_repeat) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  }
  texture.size         = size;
  texture.num_channels = num_channels;
  texture.is_srgb      = false;
  texture.is_float     = as_float;
  texture.linear       = linear;
  texture.mipmap       = mipmap;
  assert_ogl_error();
}

// cleanup
ogl_texture::~ogl_texture() { clear_texture(*this); }

// check if texture is initialized
bool is_initialized(const ogl_texture& texture) {
  return texture.texture_id != 0;
}

// clear texture
void clear_texture(ogl_texture& texture) {
  if (texture.texture_id) glDeleteTextures(1, &texture.texture_id);
  texture.size         = {0, 0};
  texture.num_channels = 0;
  texture.is_srgb      = false;
  texture.is_float     = false;
  texture.linear       = false;
  texture.mipmap       = false;
  texture.texture_id   = 0;
  assert_ogl_error();
}

void set_texture(ogl_texture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 4, (const byte*)img.data(), as_srgb,
      linear, mipmap);
}
void set_texture(ogl_texture& texture, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 4, (const float*)img.data(), as_float,
      linear, mipmap);
}

void set_texture(ogl_texture& texture, const image<vec3b>& img, bool as_srgb,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 3, (const byte*)img.data(), as_srgb,
      linear, mipmap);
}
void set_texture(ogl_texture& texture, const image<vec3f>& img, bool as_float,
    bool linear, bool mipmap) {
  set_texture(texture, img.imsize(), 3, (const float*)img.data(), as_float,
      linear, mipmap);
}

void set_cubemap(ogl_cubemap& cubemap, int size, int num_channels,
    const array<byte*, 6>& images, bool as_srgb, bool linear, bool mipmap) {
  static auto sformat = vector<uint>{
      0, GL_SRGB, GL_SRGB, GL_SRGB, GL_SRGB_ALPHA};
  static auto iformat = vector<uint>{0, GL_RGB, GL_RGB, GL_RGB, GL_RGBA};
  static auto cformat = vector<uint>{0, GL_RED, GL_RG, GL_RGB, GL_RGBA};
  assert_ogl_error();

  if (size == 0) {
    clear_cubemap(cubemap);
    return;
  }

  if (!cubemap.cubemap_id) glGenTextures(1, &cubemap.cubemap_id);
  if (cubemap.size != size || cubemap.num_channels != num_channels ||
      cubemap.is_srgb != as_srgb || cubemap.is_float == true ||
      cubemap.linear != linear || cubemap.mipmap != mipmap) {
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap.cubemap_id);

    for (auto i = 0; i < 6; i++) {
      if (!images[i]) {
        throw std::runtime_error{"cannot initialize cubemap from empty image"};
        return;
      }
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0,
          as_srgb ? sformat.at(num_channels) : iformat.at(num_channels), size,
          size, 0, cformat.at(num_channels), GL_UNSIGNED_BYTE, images[i]);
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
    // glBindTexture(GL_TEXTURE_2D, cubemap.cubemap_id);
    // glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size.x, size.y,
    //     cformat.at(num_channels), GL_UNSIGNED_BYTE, img);
    // if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
  cubemap.size         = size;
  cubemap.num_channels = num_channels;
  cubemap.is_srgb      = as_srgb;
  cubemap.is_float     = false;
  cubemap.linear       = linear;
  cubemap.mipmap       = mipmap;
  assert_ogl_error();
}

void set_cubemap(ogl_cubemap& cubemap, int size, int num_channels,
    const array<float*, 6>& images, bool as_float, bool linear, bool mipmap) {
  static auto fformat = vector<uint>{
      0, GL_RGB16F, GL_RGB16F, GL_RGB16F, GL_RGBA32F};
  static auto iformat = vector<uint>{0, GL_RGB, GL_RGB, GL_RGB, GL_RGBA};
  static auto cformat = vector<uint>{0, GL_RED, GL_RG, GL_RGB, GL_RGBA};
  assert_ogl_error();

  if (size == 0) {
    clear_cubemap(cubemap);
    return;
  }

  if (!cubemap.cubemap_id) glGenTextures(1, &cubemap.cubemap_id);
  if (cubemap.size != size || cubemap.num_channels != num_channels ||
      cubemap.is_float != as_float || cubemap.is_srgb == true ||
      cubemap.linear != linear || cubemap.mipmap != mipmap) {
    glGenTextures(1, &cubemap.cubemap_id);

    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap.cubemap_id);

    for (auto i = 0; i < 6; i++) {
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0,
          as_float ? fformat.at(num_channels) : iformat.at(num_channels), size,
          size, 0, iformat.at(num_channels), GL_FLOAT, images[i]);
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

    //    glBindTexture(GL_TEXTURE_2D, cubemap.cubemap_id);
    //    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size, size,
    //        iformat.at(num_channels), GL_FLOAT, img);
    //    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
  cubemap.size         = size;
  cubemap.num_channels = num_channels;
  cubemap.is_srgb      = false;
  cubemap.is_float     = as_float;
  cubemap.linear       = linear;
  cubemap.mipmap       = mipmap;
  assert_ogl_error();
}

// cleanup
ogl_cubemap::~ogl_cubemap() { clear_cubemap(*this); }

// check if cubemap is initialized
bool is_initialized(const ogl_cubemap& cubemap) {
  return cubemap.cubemap_id != 0;
}

// clear cubemap
void clear_cubemap(ogl_cubemap& cubemap) {
  if (cubemap.cubemap_id) glDeleteTextures(1, &cubemap.cubemap_id);
  cubemap.size         = 0;
  cubemap.num_channels = 0;
  cubemap.is_srgb      = false;
  cubemap.is_float     = false;
  cubemap.linear       = false;
  cubemap.mipmap       = false;
  cubemap.cubemap_id   = 0;
  assert_ogl_error();
}

void set_cubemap(ogl_cubemap& cubemap, const array<image<vec4b>, 6>& img,
    int num_channels, bool as_srgb, bool linear, bool mipmap) {
  auto data = array<byte*, 6>{(byte*)img[0].data(), (byte*)img[1].data(),
      (byte*)img[2].data(), (byte*)img[3].data(), (byte*)img[4].data(),
      (byte*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, num_channels, data, as_srgb, linear, mipmap);
}
void set_cubemap(ogl_cubemap& cubemap, const array<image<vec4f>, 6>& img,
    int num_channels, bool as_float, bool linear, bool mipmap) {
  auto data = array<float*, 6>{(float*)img[0].data(), (float*)img[1].data(),
      (float*)img[2].data(), (float*)img[3].data(), (float*)img[4].data(),
      (float*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, num_channels, data, as_float, linear, mipmap);
}
void set_cubemap(ogl_cubemap& cubemap, const array<image<vec3b>, 6>& img,
    int num_channels, bool as_srgb, bool linear, bool mipmap) {
  auto data = array<byte*, 6>{(byte*)img[0].data(), (byte*)img[1].data(),
      (byte*)img[2].data(), (byte*)img[3].data(), (byte*)img[4].data(),
      (byte*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, num_channels, data, as_srgb, linear, mipmap);
}
void set_cubemap(ogl_cubemap& cubemap, const array<image<vec3f>, 6>& img,
    int num_channels, bool as_float, bool linear, bool mipmap) {
  auto data = array<float*, 6>{(float*)img[0].data(), (float*)img[1].data(),
      (float*)img[2].data(), (float*)img[3].data(), (float*)img[4].data(),
      (float*)img[5].data()};
  set_cubemap(
      cubemap, img[0].imsize().x, num_channels, data, as_float, linear, mipmap);
}

// cleanup
ogl_arraybuffer::~ogl_arraybuffer() { clear_arraybuffer(*this); }

// check if buffer is initialized
bool is_initialized(const ogl_arraybuffer& buffer) {
  return buffer.buffer_id != 0;
}

// set buffer
void set_arraybuffer(ogl_arraybuffer& buffer, size_t size, int esize,
    const float* data, bool dynamic) {
  assert_ogl_error();

  if (size == 0) {
    clear_arraybuffer(buffer);
    return;
  }

  auto target = GL_ARRAY_BUFFER;
  if (size > buffer.capacity) {
    // reallocate buffer if needed
    if (buffer.buffer_id) {
      glDeleteBuffers(1, &buffer.buffer_id);
    }
    glGenBuffers(1, &buffer.buffer_id);
    glBindBuffer(target, buffer.buffer_id);
    glBufferData(target, size * sizeof(float), data,
        dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    buffer.capacity = size;
  } else {
    // we have enough space
    assert(buffer.buffer_id);
    glBindBuffer(target, buffer.buffer_id);
    glBufferSubData(target, 0, size * sizeof(float), data);
  }

  buffer.element_size = esize;
  buffer.num_elements = size / esize;
  buffer.dynamic      = dynamic;
  assert_ogl_error();
}

// clear buffer
void clear_arraybuffer(ogl_arraybuffer& buffer) {
  assert_ogl_error();
  if (buffer.buffer_id) glDeleteBuffers(1, &buffer.buffer_id);
  buffer.capacity     = 0;
  buffer.num_elements = 0;
  buffer.element_size = 0;
  buffer.dynamic      = false;
  buffer.buffer_id    = 0;
  assert_ogl_error();
}

// set buffer
void set_arraybuffer(
    ogl_arraybuffer& buffer, const vector<float>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 1, 1, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer& buffer, const vector<vec2f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 2, 2, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer& buffer, const vector<vec3f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 3, 3, (float*)data.data(), dynamic);
}
void set_arraybuffer(
    ogl_arraybuffer& buffer, const vector<vec4f>& data, bool dynamic) {
  set_arraybuffer(buffer, data.size() * 4, 4, (float*)data.data(), dynamic);
}

void set_elementbuffer(ogl_elementbuffer& buffer, size_t size, int esize,
    const int* data, bool dynamic) {
  assert_ogl_error();

  if (size == 0) {
    clear_elementbuffer(buffer);
    return;
  }

  auto target = GL_ELEMENT_ARRAY_BUFFER;
  if (size > buffer.capacity) {
    // reallocate buffer if needed
    if (buffer.buffer_id) {
      glDeleteBuffers(1, &buffer.buffer_id);
    }
    glGenBuffers(1, &buffer.buffer_id);
    glBindBuffer(target, buffer.buffer_id);
    glBufferData(target, size * sizeof(int), data,
        dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    buffer.capacity = size;
  } else {
    // we have enough space
    assert(buffer.buffer_id);
    glBindBuffer(target, buffer.buffer_id);
    glBufferSubData(target, 0, size * sizeof(int), data);
  }

  buffer.element_size = esize;
  buffer.num_elements = size / esize;
  buffer.dynamic      = dynamic;
  assert_ogl_error();
}

// cleanup
ogl_elementbuffer::~ogl_elementbuffer() { clear_elementbuffer(*this); }

// check if buffer is initialized
bool is_initialized(const ogl_elementbuffer& buffer) {
  return buffer.buffer_id != 0;
}

// clear buffer
void clear_elementbuffer(ogl_elementbuffer& buffer) {
  assert_ogl_error();
  if (buffer.buffer_id) glDeleteBuffers(1, &buffer.buffer_id);
  buffer.capacity     = 0;
  buffer.num_elements = 0;
  buffer.element_size = 0;
  buffer.dynamic      = false;
  buffer.buffer_id    = 0;
  assert_ogl_error();
}

// set buffer
void set_elementbuffer(
    ogl_elementbuffer& buffer, const vector<int>& points, bool dynamic) {
  set_elementbuffer(buffer, points.size() * 1, 1, (int*)points.data(), dynamic);
}
void set_elementbuffer(
    ogl_elementbuffer& buffer, const vector<vec2i>& lines, bool dynamic) {
  set_elementbuffer(buffer, lines.size() * 2, 2, (int*)lines.data(), dynamic);
}
void set_elementbuffer(
    ogl_elementbuffer& buffer, const vector<vec3i>& triangles, bool dynamic) {
  set_elementbuffer(
      buffer, triangles.size() * 3, 3, (int*)triangles.data(), dynamic);
}

// initialize program
bool set_program(ogl_program& program, const string& vertex,
    const string& fragment, string& error, string& errorlog) {
  // error
  auto program_error = [&error, &errorlog, &program](
                           const char* message, const char* log) {
    clear_program(program);
    error    = message;
    errorlog = log;
    return false;
  };

  // clear
  if (program.program_id) clear_program(program);

  // setup code
  program.vertex_code   = vertex;
  program.fragment_code = fragment;

  const char* ccvertex   = vertex.data();
  const char* ccfragment = fragment.data();
  auto        errflags   = 0;
  auto        errbuf     = array<char, 10000>{};

  // create vertex
  assert_ogl_error();
  program.vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(program.vertex_id, 1, &ccvertex, NULL);
  glCompileShader(program.vertex_id);
  glGetShaderiv(program.vertex_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(program.vertex_id, 10000, 0, errbuf.data());
    return program_error("vertex shader not compiled", errbuf.data());
  }
  assert_ogl_error();

  // create fragment
  assert_ogl_error();
  program.fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(program.fragment_id, 1, &ccfragment, NULL);
  glCompileShader(program.fragment_id);
  glGetShaderiv(program.fragment_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(program.fragment_id, 10000, 0, errbuf.data());
    return program_error("fragment shader not compiled", errbuf.data());
  }
  assert_ogl_error();

  // create program
  assert_ogl_error();
  program.program_id = glCreateProgram();
  glAttachShader(program.program_id, program.vertex_id);
  glAttachShader(program.program_id, program.fragment_id);
  glLinkProgram(program.program_id);
  glGetProgramiv(program.program_id, GL_LINK_STATUS, &errflags);
  if (errflags == 0) {
    glGetProgramInfoLog(program.program_id, 10000, 0, errbuf.data());
    return program_error("program not linked", errbuf.data());
  }
  // TODO(giacomo): Apparently validation must be done just before drawing.
  //    https://community.khronos.org/t/samplers-of-different-types-use-the-same-textur/66329
  // If done here, validation fails when using cubemaps and textures in the
  // same shader. We should create a function validate_program() anc call it
  // separately.
  //
  // glValidateProgram(program.program_id);
  // glGetProgramiv(program.program_id, GL_VALIDATE_STATUS, &errflags);
  // if (!errflags) {
  //   glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
  //   return program_error("program not validated", errbuf);
  // }
  assert_ogl_error();

  // done
  return true;
}

// initialize program
bool set_program(ogl_program& program, const string& vertex,
    const string& fragment, string& error) {
  auto errorlog = string{};
  return set_program(program, vertex, fragment, error, errorlog);
}

// initialize program, print eventual errors to stdout
bool set_program(ogl_program& program, const string& vertex,
    const string& fragment, bool exceptions) {
  auto error    = string{};
  auto errorlog = string{};
  if (!set_program(program, vertex, fragment, error, errorlog)) {
    if (exceptions) throw std::runtime_error{error + "\n" + errorlog};
    return false;
  }
  return true;
}

// clear program
void clear_program(ogl_program& program) {
  if (program.program_id) glDeleteProgram(program.program_id);
  if (program.vertex_id) glDeleteShader(program.vertex_id);
  if (program.fragment_id) glDeleteShader(program.fragment_id);
  program.vertex_code   = {};
  program.fragment_code = {};
  program.program_id    = 0;
  program.vertex_id     = 0;
  program.fragment_id   = 0;
  assert_ogl_error();
}

// cleanup
ogl_program::~ogl_program() { clear_program(*this); }

bool is_initialized(const ogl_program& program) {
  return program.program_id != 0;
}

// bind program
void bind_program(const ogl_program& program) {
  assert_ogl_error();
  glUseProgram(program.program_id);
  ogl_program::bound_program_id = program.program_id;
  assert_ogl_error();
}

// unbind program
void unbind_program() {
  glUseProgram(0);
  assert_ogl_error();
  ogl_program::bound_program_id = 0;
}

// set uniforms
void set_uniform(const ogl_program& program, int location, int value) {
  glUniform1i(location, value);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec2i& value) {
  glUniform2i(location, value.x, value.y);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec3i& value) {
  glUniform3i(location, value.x, value.y, value.z);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec4i& value) {
  glUniform4i(location, value.x, value.y, value.z, value.w);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, float value) {
  glUniform1f(location, value);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec2f& value) {
  glUniform2f(location, value.x, value.y);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec3f& value) {
  glUniform3f(location, value.x, value.y, value.z);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const vec4f& value) {
  glUniform4f(location, value.x, value.y, value.z, value.w);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const mat2f& value) {
  glUniformMatrix2fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const mat3f& value) {
  glUniformMatrix3fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(const ogl_program& program, int location, const mat4f& value) {
  glUniformMatrix4fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(
    const ogl_program& program, int location, const frame2f& value) {
  glUniformMatrix3x2fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

void set_uniform(
    const ogl_program& program, int location, const frame3f& value) {
  glUniformMatrix4x3fv(location, 1, false, &value.x.x);
  assert_ogl_error();
}

// get uniform location
int get_uniform_location(const ogl_program& program, const char* name) {
  return glGetUniformLocation(program.program_id, name);
}

// set uniform texture
void set_uniform(const ogl_program& program, int location,
    const ogl_texture& texture, int unit) {
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glUniform1i(location, unit);
  assert_ogl_error();
}
void set_uniform(const ogl_program& program, const char* name,
    const ogl_texture& texture, int unit) {
  return set_uniform(
      program, get_uniform_location(program, name), texture, unit);
}
void set_uniform(const ogl_program& program, int location, int location_on,
    const ogl_texture& texture, int unit) {
  assert_ogl_error();
  if (texture.texture_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glUniform1i(location, unit);
    glUniform1i(location_on, 1);
  } else {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, 0);
    glUniform1i(location, unit);
    glUniform1i(location_on, 0);
  }
  assert_ogl_error();
}
void set_uniform(const ogl_program& program, const char* name,
    const char* name_on, const ogl_texture& texture, int unit) {
  return set_uniform(program, get_uniform_location(program, name),
      get_uniform_location(program, name_on), texture, unit);
}

// set uniform cubemap
void set_uniform(const ogl_program& program, int location,
    const ogl_cubemap& cubemap, int unit) {
  assert_ogl_error();
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap.cubemap_id);
  glUniform1i(location, unit);
  assert_ogl_error();
}
void set_uniform(const ogl_program& program, const char* name,
    const ogl_cubemap& cubemap, int unit) {
  return set_uniform(
      program, get_uniform_location(program, name), cubemap, unit);
}
void set_uniform(const ogl_program& program, int location, int location_on,
    const ogl_cubemap& cubemap, int unit) {
  assert_ogl_error();
  if (cubemap.cubemap_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap.cubemap_id);
    glUniform1i(location, unit);
    glUniform1i(location_on, 1);
  } else {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
    glUniform1i(location, unit);
    glUniform1i(location_on, 0);
  }
  assert_ogl_error();
}
void set_uniform(const ogl_program& program, const char* name,
    const char* name_on, const ogl_cubemap& cubemap, int unit) {
  return set_uniform(program, get_uniform_location(program, name),
      get_uniform_location(program, name_on), cubemap, unit);
}

// cleanup
ogl_framebuffer::~ogl_framebuffer() { clear_framebuffer(*this); }

void set_framebuffer(ogl_framebuffer& framebuffer, const vec2i& size) {
  if (size == zero2i) {
    clear_framebuffer(framebuffer);
    return;
  }

  if (!framebuffer.framebuffer_id) {
    glGenFramebuffers(1, &framebuffer.framebuffer_id);
  }

  if (!framebuffer.renderbuffer_id) {
    glGenRenderbuffers(1, &framebuffer.renderbuffer_id);
    // bind together frame buffer and render buffer
    // TODO(giacomo): Why do we need to put STENCIL8 to make things work on
    // Mac??
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.framebuffer_id);
    glBindRenderbuffer(GL_RENDERBUFFER, framebuffer.renderbuffer_id);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT,
        GL_RENDERBUFFER, framebuffer.renderbuffer_id);
  }

  if (size != framebuffer.size) {
    // create render buffer for depth and stencil
    // TODO(giacomo): We put STENCIL here for the same reason...
    glBindRenderbuffer(GL_RENDERBUFFER, framebuffer.renderbuffer_id);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, size.x, size.y);
    framebuffer.size = size;
  }

  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.framebuffer_id);
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
  glBindFramebuffer(GL_FRAMEBUFFER, ogl_framebuffer::bound_framebuffer_id);

  assert_ogl_error();
}

inline void set_framebuffer_texture(const ogl_framebuffer& framebuffer,
    uint texture_id, uint target, uint mipmap_level) {
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.framebuffer_id);
  glFramebufferTexture2D(
      GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, target, texture_id, mipmap_level);
  glBindFramebuffer(GL_FRAMEBUFFER, ogl_framebuffer::bound_framebuffer_id);
  assert_ogl_error();
}

bool is_framebuffer_bound(const ogl_framebuffer& framebuffer) {
  return framebuffer.framebuffer_id == ogl_framebuffer::bound_framebuffer_id;
}

void set_framebuffer_texture(const ogl_framebuffer& framebuffer,
    const ogl_texture& texture, uint mipmap_level) {
  set_framebuffer_texture(
      framebuffer, texture.texture_id, GL_TEXTURE_2D, mipmap_level);
}

void set_framebuffer_texture(const ogl_framebuffer& framebuffer,
    const ogl_cubemap& cubemap, uint face, uint mipmap_level) {
  set_framebuffer_texture(framebuffer, cubemap.cubemap_id,
      GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, mipmap_level);
}

void bind_framebuffer(const ogl_framebuffer& framebuffer) {
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.framebuffer_id);
  assert_ogl_error();
  ogl_framebuffer::bound_framebuffer_id = framebuffer.framebuffer_id;
}

void unbind_framebuffer() {
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  assert_ogl_error();
  ogl_framebuffer::bound_framebuffer_id = 0;
}

void clear_framebuffer(ogl_framebuffer& framebuffer) {
  if (ogl_framebuffer::bound_framebuffer_id == framebuffer.framebuffer_id) {
    unbind_framebuffer();
  }
  glDeleteFramebuffers(1, &framebuffer.framebuffer_id);
  glDeleteRenderbuffers(1, &framebuffer.renderbuffer_id);
  framebuffer.size            = {0, 0};
  framebuffer.framebuffer_id  = 0;
  framebuffer.renderbuffer_id = 0;
  assert_ogl_error();
}

ogl_shape::~ogl_shape() { clear_shape(*this); }

void bind_shape(const ogl_shape& shape) { glBindVertexArray(shape.shape_id); }

bool is_initialized(const ogl_shape& shape) { return shape.shape_id != 0; }

// Clear an OpenGL shape
void clear_shape(ogl_shape& shape) {
  for (auto& buffer : shape.vertex_buffers) {
    clear_arraybuffer(buffer);
  }
  clear_elementbuffer(shape.index_buffer);
  glDeleteVertexArrays(1, &shape.shape_id);
  shape.num_instances = 0;
  shape.shape_id      = 0;
  assert_ogl_error();
}

template <typename T>
void set_vertex_buffer_impl(
    ogl_shape& shape, const vector<T>& data, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  while (shape.vertex_buffers.size() <= location) {
    shape.vertex_buffers.emplace_back();
  }
  set_arraybuffer(shape.vertex_buffers[location], data, false);
  glBindVertexArray(shape.shape_id);
  auto& buffer = shape.vertex_buffers[location];
  assert_ogl_error();
  glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
  glEnableVertexAttribArray(location);
  glVertexAttribPointer(
      location, buffer.element_size, GL_FLOAT, false, 0, nullptr);
  assert_ogl_error();
}

void set_vertex_buffer(
    ogl_shape& shape, const vector<float>& values, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  set_vertex_buffer_impl(shape, values, location);
}
void set_vertex_buffer(
    ogl_shape& shape, const vector<vec2f>& values, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  set_vertex_buffer_impl(shape, values, location);
}
void set_vertex_buffer(
    ogl_shape& shape, const vector<vec3f>& values, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  set_vertex_buffer_impl(shape, values, location);
}
void set_vertex_buffer(
    ogl_shape& shape, const vector<vec4f>& values, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  set_vertex_buffer_impl(shape, values, location);
}

void set_vertex_buffer(ogl_shape& shape, float value, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  glBindVertexArray(shape.shape_id);
  glVertexAttrib1f(location, value);
  assert_ogl_error();
}
void set_vertex_buffer(ogl_shape& shape, const vec2f& value, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  glBindVertexArray(shape.shape_id);
  glVertexAttrib2f(location, value.x, value.y);
  assert_ogl_error();
}
void set_vertex_buffer(ogl_shape& shape, const vec3f& value, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  glBindVertexArray(shape.shape_id);
  glVertexAttrib3f(location, value.x, value.y, value.z);
  assert_ogl_error();
}
void set_vertex_buffer(ogl_shape& shape, const vec4f& value, int location) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  glBindVertexArray(shape.shape_id);
  glVertexAttrib4f(location, value.x, value.y, value.z, value.w);
  assert_ogl_error();
}

void set_instance_buffer(ogl_shape& shape, int location, bool is_instance) {
  if (!shape.shape_id) glGenVertexArrays(1, &shape.shape_id);
  bind_shape(shape);
  if (is_instance) {
    glVertexAttribDivisor(location, 1);
    shape.num_instances = shape.vertex_buffers[location].num_elements;
  } else {
    glVertexAttribDivisor(location, 0);
    shape.num_instances = 0;
  }
  assert_ogl_error();
}

void set_index_buffer(ogl_shape& shape, const vector<int>& indices) {
  set_elementbuffer(shape.index_buffer, indices);
  shape.elements = ogl_element_type::points;
}
void set_index_buffer(ogl_shape& shape, const vector<vec2i>& indices) {
  set_elementbuffer(shape.index_buffer, indices);
  shape.elements = ogl_element_type::lines;
}
void set_index_buffer(ogl_shape& shape, const vector<vec3i>& indices) {
  set_elementbuffer(shape.index_buffer, indices);
  shape.elements = ogl_element_type::triangles;
}

// set point size
void set_point_size(ogl_shape& shape, float point_size) {
  shape.point_size = point_size;
}

void draw_shape(const ogl_shape& shape) {
  if (shape.shape_id == 0) return;
  bind_shape(shape);
  auto type = GL_TRIANGLES;
  switch (shape.elements) {
    case ogl_element_type::points: type = GL_POINTS; break;
    case ogl_element_type::lines: type = GL_LINES; break;
    case ogl_element_type::line_strip: type = GL_LINE_STRIP; break;
    case ogl_element_type::triangles: type = GL_TRIANGLES; break;
    case ogl_element_type::triangle_strip: type = GL_TRIANGLE_STRIP; break;
    case ogl_element_type::triangle_fan: type = GL_TRIANGLE_FAN; break;
  }

  if (shape.elements == ogl_element_type::points) {
    glPointSize(shape.point_size);
  }

  auto& indices = shape.index_buffer;
  if (indices.buffer_id != 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices.buffer_id);
    if (shape.num_instances == 0) {
      glDrawElements(type, (GLsizei)indices.num_elements * indices.element_size,
          GL_UNSIGNED_INT, nullptr);
    } else {
      glDrawElementsInstanced(type,
          (GLsizei)indices.num_elements * indices.element_size, GL_UNSIGNED_INT,
          nullptr, (GLsizei)shape.num_instances);
    }
  } else {
    auto& vertices = shape.vertex_buffers[0];
    glDrawArrays(type, 0, (int)vertices.num_elements);
  }

  if (shape.elements == ogl_element_type::points) {
    glPointSize(shape.point_size);
  }

  assert_ogl_error();
}

void set_cube_shape(ogl_shape& shape) {
  // clang-format off
  static const auto positions = vector<vec3f>{
    {1, -1, -1}, {1, -1,  1}, {-1, -1,  1}, {-1, -1, -1},
    {1,  1, -1}, {1,  1,  1}, {-1,  1,  1}, {-1,  1, -1},
  };
  static const auto triangles = vector<vec3i>{
    {1, 3, 0}, {7, 5, 4}, {4, 1, 0}, {5, 2, 1},
    {2, 7, 3}, {0, 7, 4}, {1, 2, 3}, {7, 6, 5},
    {4, 5, 1}, {5, 6, 2}, {2, 6, 7}, {0, 3, 7}
  };
  // clang-format on
  set_vertex_buffer(shape, positions, 0);
  set_index_buffer(shape, triangles);
}

void set_quad_shape(ogl_shape& shape) {
  // clang-format off
  static const auto positions = vector<vec3f>{
    {-1, -1, 0}, {1, -1,  0}, {1, 1,  0}, {-1, 1, 0},
  };
  static const auto triangles = vector<vec3i>{
    {0, 1, 3}, {3, 2, 1}
  };
  // clang-format on
  set_vertex_buffer(shape, positions, 0);
  set_index_buffer(shape, triangles);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

static auto ogl_image_vertex =
    R"(
#version 330
in vec2 positions;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * image_size * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
    frag_texcoord = positions * 0.5 + 0.5;
}
)";
#if 0
static auto ogl_image_vertex = R"(
#version 330
in vec2 positions;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size, border_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * (image_size + border_size*2) * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
    frag_texcoord = positions * 0.5 + 0.5;
}
)";
#endif
static auto ogl_image_fragment =
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
static auto ogl_image_fragment = R"(
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

bool is_initialized(const ogl_image& oimg) {
  return is_initialized(oimg.program);
}

// init image program
bool init_image(ogl_image& oimg) {
  if (is_initialized(oimg)) return true;
  if (!set_program(oimg.program, ogl_image_vertex, ogl_image_fragment))
    return false;
  set_quad_shape(oimg.quad);
  return true;
}

// clear an opengl image
void clear_image(ogl_image& oimg) {
  clear_program(oimg.program);
  clear_texture(oimg.texture);
  clear_shape(oimg.quad);
}

// update image data
void set_image(
    ogl_image& oimg, const image<vec4f>& img, bool linear, bool mipmap) {
  set_texture(oimg.texture, img, false, linear, mipmap);
}
void set_image(
    ogl_image& oimg, const image<vec4b>& img, bool linear, bool mipmap) {
  set_texture(oimg.texture, img, false, linear, mipmap);
}

void set_image(
    ogl_image& oimg, const image_data& img, bool linear, bool mipmap) {
  if (is_float(img)) {
    set_texture(oimg.texture, {img.width, img.height}, 4,
        (const float*)img.pixelsf.data(), false, linear, mipmap);
  } else {
    set_texture(oimg.texture, {img.width, img.height}, 4,
        (const byte*)img.pixelsb.data(), false, linear, mipmap);
  }
}

// draw image
void draw_image(ogl_image& oimg, const ogl_image_params& params) {
  assert_ogl_error();
  set_ogl_viewport(params.framebuffer);
  clear_ogl_framebuffer(params.background);
  bind_program(oimg.program);
  set_uniform(oimg.program, "txt", oimg.texture, 0);
  set_uniform(oimg.program, "window_size",
      vec2f{(float)params.window.x, (float)params.window.y});
  set_uniform(oimg.program, "image_size",
      vec2f{(float)oimg.texture.size.x, (float)oimg.texture.size.y});
  set_uniform(oimg.program, "image_center", params.center);
  set_uniform(oimg.program, "image_scale", params.scale);
  draw_shape(oimg.quad);
  unbind_program();
}

}  // namespace yocto
