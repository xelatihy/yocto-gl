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
#include <stdarg.h>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"

#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"

#define CUTE_FILES_IMPLEMENTATION
#include "ext/cute_files.h"

namespace yocto {

void check_glerror() {
  if (glGetError() != GL_NO_ERROR) printf("gl error\n");
}

void clear_glframebuffer(const vec4f& color, bool clear_depth) {
  glClearColor(color.x, color.y, color.z, color.w);
  if (clear_depth) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
  } else {
    glClear(GL_COLOR_BUFFER_BIT);
  }
}

void set_glviewport(const vec4i& viewport) {
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);
}

void set_glwireframe(bool enabled) {
  if (enabled)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void set_glblending(bool enabled) {
  if (enabled) {
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
  } else {
    glDisable(GL_BLEND);
  }
}

void init_glprogram(
    opengl_program& program, const char* vertex, const char* fragment) {
  assert(glGetError() == GL_NO_ERROR);
  glGenVertexArrays(1, &program.vertex_array_object_id);
  glBindVertexArray(program.vertex_array_object_id);
  assert(glGetError() == GL_NO_ERROR);

  int  errflags;
  char errbuf[10000];

  // create vertex
  assert(glGetError() == GL_NO_ERROR);
  program.vertex_shader_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(program.vertex_shader_id, 1, &vertex, NULL);
  glCompileShader(program.vertex_shader_id);
  glGetShaderiv(program.vertex_shader_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program.vertex_shader_id, 10000, 0, errbuf);
    throw gl_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create fragment
  assert(glGetError() == GL_NO_ERROR);
  program.fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(program.fragment_shader_id, 1, &fragment, NULL);
  glCompileShader(program.fragment_shader_id);
  glGetShaderiv(program.fragment_shader_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program.fragment_shader_id, 10000, 0, errbuf);
    throw gl_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create program
  assert(glGetError() == GL_NO_ERROR);
  program.program_id = glCreateProgram();
  glAttachShader(program.program_id, program.vertex_shader_id);
  glAttachShader(program.program_id, program.fragment_shader_id);
  glLinkProgram(program.program_id);
  glValidateProgram(program.program_id);
  glGetProgramiv(program.program_id, GL_LINK_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
    throw gl_error("program not linked with error\n"s + errbuf);
  }
  glGetProgramiv(program.program_id, GL_VALIDATE_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
    throw gl_error("program not linked with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void delete_glprogram(opengl_program& program) {
  if (!program) return;
  glDeleteProgram(program.program_id);
  glDeleteShader(program.vertex_shader_id);
  glDeleteShader(program.fragment_shader_id);
  program = {};
}

void init_gltexture(opengl_texture& texture, const vec2i& size, bool as_float,
    bool as_srgb, bool linear, bool mipmap) {
  texture = opengl_texture();
  assert(glGetError() == GL_NO_ERROR);
  glGenTextures(1, &texture.texture_id);
  texture.size = size;
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  if (as_float) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  } else if (as_srgb) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, size.x, size.y, 0, GL_RGBA,
        GL_UNSIGNED_BYTE, nullptr);
  } else {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  }
  if (mipmap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture(
    opengl_texture& texture, const image<vec4f>& img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
      GL_FLOAT, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture_region(opengl_texture& texture, const image<vec4f>& img,
    const image_region& region, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  auto clipped = image<vec4f>{};
  get_region(clipped, img, region);
  glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y, region.size().x,
      region.size().y, GL_RGBA, GL_FLOAT, clipped.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture(
    opengl_texture& texture, const image<vec4b>& img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
      GL_UNSIGNED_BYTE, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture_region(opengl_texture& texture, const image<vec4b>& img,
    const image_region& region, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  auto clipped = image<vec4b>{};
  get_region(clipped, img, region);
  glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y, region.size().x,
      region.size().y, GL_RGBA, GL_UNSIGNED_BYTE, clipped.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void delete_gltexture(opengl_texture& texture) {
  if (!texture) return;
  glDeleteTextures(1, &texture.texture_id);
  texture = {};
}

template <typename T>
void init_glarray_buffer_impl(
    opengl_arraybuffer& buffer, const vector<T>& array, bool dynamic) {
  buffer           = opengl_arraybuffer{};
  buffer.num       = size(array);
  buffer.elem_size = sizeof(T);
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer.buffer_id);
  glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
  glBufferData(GL_ARRAY_BUFFER, size(array) * sizeof(T), array.data(),
      (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<float>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec2f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec3f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec4f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}

void delete_glarraybuffer(opengl_arraybuffer& buffer) {
  if (!buffer) return;
  glDeleteBuffers(1, &buffer.buffer_id);
  buffer = {};
}

template <typename T>
void init_glelementbuffer_impl(
    opengl_elementbuffer& buffer, const vector<T>& array, bool dynamic) {
  buffer           = opengl_elementbuffer{};
  buffer.num       = size(array);
  buffer.elem_size = sizeof(T);
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer.buffer_id);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, size(array) * sizeof(T), array.data(),
      (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<int>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}
void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<vec2i>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}
void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<vec3i>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}

void delete_glelementbuffer(opengl_elementbuffer& buffer) {
  if (!buffer) return;
  glDeleteBuffers(1, &buffer.buffer_id);
  buffer = {};
}

void bind_glprogram(opengl_program& program) {
  glUseProgram(program.program_id);
}
void unbind_opengl_program() { glUseProgram(0); }

int get_gluniform_location(const opengl_program& program, const char* name) {
  return glGetUniformLocation(program.program_id, name);
}

void set_gluniform(int locatiom, int value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform1i(locatiom, value);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec2i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform2i(locatiom, value.x, value.y);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec3i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform3i(locatiom, value.x, value.y, value.z);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec4i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform4i(locatiom, value.x, value.y, value.z, value.w);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, float value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform1f(locatiom, value);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec2f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform2f(locatiom, value.x, value.y);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform3f(locatiom, value.x, value.y, value.z);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform4f(locatiom, value.x, value.y, value.z, value.w);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const mat4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniformMatrix4fv(locatiom, 1, false, &value.x.x);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const frame3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniformMatrix4x3fv(locatiom, 1, false, &value.x.x);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(
    int locatiom, const opengl_texture& texture, int unit) {
  assert(glGetError() == GL_NO_ERROR);
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glUniform1i(locatiom, unit);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(opengl_program& program, const char* name,
    const opengl_texture& texture, int unit) {
  set_gluniform_texture(get_gluniform_location(program, name), texture, unit);
}

void set_gluniform_texture(
    int locatiom, int locatiom_on, const opengl_texture& texture, int unit) {
  assert(glGetError() == GL_NO_ERROR);
  if (texture.texture_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glUniform1i(locatiom, unit);
    glUniform1i(locatiom_on, 1);
  } else {
    glUniform1i(locatiom_on, 0);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(opengl_program& program, const char* name,
    const char* name_on, const opengl_texture& texture, int unit) {
  set_gluniform_texture(get_gluniform_location(program, name),
      get_gluniform_location(program, name_on), texture, unit);
}

int get_glvertexattrib_location(
    const opengl_program& program, const char* name) {
  return glGetAttribLocation(program.program_id, name);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, float value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 1, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib1f(locatiom, value);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec2f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 2, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib2f(locatiom, value.x, value.y);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 3, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib3f(locatiom, value.x, value.y, value.z);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 4, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib4f(locatiom, value.x, value.y, value.z, value.w);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void draw_glpoints(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_POINTS, num, GL_UNSIGNED_INT, nullptr);
}

void draw_gllines(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_LINES, num * 2, GL_UNSIGNED_INT, nullptr);
}

void draw_gltriangles(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

void draw_glimage(const opengl_texture& texture, int win_width, int win_height,
    const vec2f& image_center, float image_scale) {
  static opengl_program       gl_prog      = {};
  static opengl_arraybuffer   gl_texcoord  = {};
  static opengl_elementbuffer gl_triangles = {};

  // initialization
  if (!gl_prog) {
    auto vert = R"(
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
    auto frag = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform sampler2D txt;
            void main() {
                frag_color = texture(txt, frag_texcoord);
            }
        )";
    init_glprogram(gl_prog, vert, frag);
    init_glarraybuffer(
        gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    init_glelementbuffer(
        gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
  }

  // draw
  check_glerror();
  bind_glprogram(gl_prog);
  set_gluniform_texture(gl_prog, "txt", texture, 0);
  set_gluniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_gluniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_gluniform(gl_prog, "image_center", image_center);
  set_gluniform(gl_prog, "image_scale", image_scale);
  set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
  draw_gltriangles(gl_triangles, 2);
  unbind_opengl_program();
  check_glerror();
}

void draw_glimage_background(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size) {
  static opengl_program       gl_prog      = {};
  static opengl_arraybuffer   gl_texcoord  = {};
  static opengl_elementbuffer gl_triangles = {};

  // initialization
  if (!gl_prog) {
    auto vert = R"(
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
    auto frag = R"(
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
    init_glprogram(gl_prog, vert, frag);
    init_glarraybuffer(
        gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    init_glelementbuffer(
        gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
  }

  // draw
  bind_glprogram(gl_prog);
  set_gluniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_gluniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_gluniform(
      gl_prog, "border_size", vec2f{(float)border_size, (float)border_size});
  set_gluniform(gl_prog, "image_center", image_center);
  set_gluniform(gl_prog, "image_scale", image_scale);
  set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
  draw_gltriangles(gl_triangles, 2);
  unbind_opengl_program();
}

void _glfw_refresh_callback(GLFWwindow* glfw) {
  auto& win = *(const opengl_window*)glfwGetWindowUserPointer(glfw);
  if (win.refresh_cb) win.refresh_cb(win);
}

void _glfw_drop_callback(GLFWwindow* glfw, int num, const char** paths) {
  auto& win = *(const opengl_window*)glfwGetWindowUserPointer(glfw);
  if (win.drop_cb) {
    auto pathv = vector<string>();
    for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
    win.drop_cb(win, pathv);
  }
}

void init_glwindow(opengl_window& win, const vec2i& size, const string& title,
    void* user_pointer, refresh_glcallback refresh_cb) {
  // init glfw
  if (!glfwInit()) throw gl_error("cannot initialize windowing system");
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // create window
  win     = opengl_window();
  win.win = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
  if (!win.win) throw gl_error("cannot initialize windowing system");
  glfwMakeContextCurrent(win.win);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowRefreshCallback(win.win, _glfw_refresh_callback);
  glfwSetWindowUserPointer(win.win, &win);
  win.user_ptr   = user_pointer;
  win.refresh_cb = refresh_cb;

  // init gl extensions
  if (!gladLoadGL()) throw gl_error("cannot initialize OpenGL extensions");
}

void delete_glwindow(opengl_window& win) {
  glfwDestroyWindow(win.win);
  glfwTerminate();
  win.win = nullptr;
}

void* get_gluser_pointer(const opengl_window& win) { return win.user_ptr; }

void set_drop_glcallback(opengl_window& win, drop_glcallback drop_cb) {
  win.drop_cb = drop_cb;
  glfwSetDropCallback(win.win, _glfw_drop_callback);
}

vec2i get_glframebuffer_size(const opengl_window& win, bool ignore_widgets) {
  auto size = zero2i;
  glfwGetFramebufferSize(win.win, &size.x, &size.y);
  if (ignore_widgets && win.widgets_width) {
    auto win_size = zero2i;
    glfwGetWindowSize(win.win, &win_size.x, &win_size.y);
    size.x -= (int)(win.widgets_width * (float)size.x / (float)win_size.x);
  }
  return size;
}

vec4i get_glframebuffer_viewport(
    const opengl_window& win, bool ignore_widgets) {
  auto viewport = zero4i;
  glfwGetFramebufferSize(win.win, &viewport.z, &viewport.w);
  if (ignore_widgets && win.widgets_width) {
    auto win_size = zero2i;
    glfwGetWindowSize(win.win, &win_size.x, &win_size.y);
    auto offset = (int)(win.widgets_width * (float)viewport.z / win_size.x);
    viewport.z -= offset;
    if (win.widgets_left) viewport.x += offset;
  }
  return viewport;
}

vec2i get_glwindow_size(const opengl_window& win, bool ignore_widgets) {
  auto size = zero2i;
  glfwGetWindowSize(win.win, &size.x, &size.y);
  if (ignore_widgets && win.widgets_width) size.x -= win.widgets_width;
  return size;
}

bool should_glwindow_close(const opengl_window& win) {
  return glfwWindowShouldClose(win.win);
}
void set_glwindow_close(const opengl_window& win, bool close) {
  glfwSetWindowShouldClose(win.win, close ? GLFW_TRUE : GLFW_FALSE);
}

vec2f get_glmouse_pos(const opengl_window& win, bool ignore_widgets) {
  double mouse_posx, mouse_posy;
  glfwGetCursorPos(win.win, &mouse_posx, &mouse_posy);
  auto pos = vec2f{(float)mouse_posx, (float)mouse_posy};
  if (ignore_widgets && win.widgets_width && win.widgets_left) {
    pos.x -= win.widgets_width;
  }
  return pos;
}

bool get_glmouse_left(const opengl_window& win) {
  return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
}
bool get_glmouse_right(const opengl_window& win) {
  return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
}

bool get_glalt_key(const opengl_window& win) {
  return glfwGetKey(win.win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
         glfwGetKey(win.win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

bool get_glshift_key(const opengl_window& win) {
  return glfwGetKey(win.win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
         glfwGetKey(win.win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

void process_glevents(const opengl_window& win, bool wait) {
  if (wait)
    glfwWaitEvents();
  else
    glfwPollEvents();
}

void swap_glbuffers(const opengl_window& win) { glfwSwapBuffers(win.win); }

void init_glwidgets(opengl_window& win, int width, bool left) {
  // init widgets
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(win.win, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();
  win.widgets_width = width;
  win.widgets_left  = left;
}

bool get_glwidgets_active(const opengl_window& win) {
  auto io = &ImGui::GetIO();
  return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

void begin_glwidgets(const opengl_window& win) {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  auto win_size = get_glwindow_size(win, false);
  if (win.widgets_left) {
    ImGui::SetNextWindowPos({0, 0});
    ImGui::SetNextWindowSize({(float)win.widgets_width, (float)win_size.y});
  } else {
    ImGui::SetNextWindowPos({(float)(win_size.x - win.widgets_width), 0});
    ImGui::SetNextWindowSize({(float)win.widgets_width, (float)win_size.y});
  }
  ImGui::SetNextWindowCollapsed(false);
  ImGui::SetNextWindowBgAlpha(1);
}

void end_glwidgets(const opengl_window& win) {
  ImGui::End();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool begin_glwidgets_window(const opengl_window& win, const char* title) {
  return ImGui::Begin(title, nullptr,
      // ImGuiWindowFlags_NoTitleBar |
      ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
          ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings);
}

bool begin_glheader(const opengl_window& win, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_glheader(const opengl_window& win) { ImGui::PopID(); }

bool begin_gltabbar(const opengl_window& win, const char* lbl) {
  if (!ImGui::BeginTabBar(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_gltabbar(const opengl_window& win) {
  ImGui::PopID();
  ImGui::EndTabBar();
}

bool begin_gltabitem(const opengl_window& win, const char* lbl) {
  if (!ImGui::BeginTabItem(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_gltabitem(const opengl_window& win) {
  ImGui::PopID();
  ImGui::EndTabItem();
}

void open_glmodal(const opengl_window& win, const char* lbl) {
  ImGui::OpenPopup(lbl);
}
void clear_glmodal(const opengl_window& win) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(const opengl_window& win, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(const opengl_window& win) { ImGui::EndPopup(); }
bool is_glmodal_open(const opengl_window& win, const char* lbl) {
  return ImGui::IsPopupOpen(lbl);
}

bool draw_glmessage(
    const opengl_window& win, const char* lbl, const string& message) {
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
bool draw_glfiledialog(const opengl_window& win, const char* lbl, string& path,
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

bool draw_glbutton(const opengl_window& win, const char* lbl) {
  return ImGui::Button(lbl);
}
bool draw_glbutton(const opengl_window& win, const char* lbl, bool enabled) {
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

void draw_gltext(const opengl_window& win, const string& text) {
  ImGui::Text("%s", text.c_str());
}

void draw_gllabel(
    const opengl_window& win, const char* lbl, const string& texture) {
  ImGui::LabelText(lbl, "%s", texture.c_str());
}

void draw_gllabel(
    const opengl_window& win, const char* lbl, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  ImGui::LabelTextV(lbl, fmt, args);
  va_end(args);
}

void draw_glseparator(const opengl_window& win) { ImGui::Separator(); }

void continue_glline(const opengl_window& win) { ImGui::SameLine(); }

bool draw_gltextinput(
    const opengl_window& win, const char* lbl, string& value) {
  char buffer[4096];
  auto num = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer, sizeof(buffer));
  if (edited) value = buffer;
  return edited;
}

bool draw_glslider(const opengl_window& win, const char* lbl, float& value,
    float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec2f& value,
    float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec3f& value,
    float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec4f& value,
    float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_glslider(
    const opengl_window& win, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_gldragger(const opengl_window& win, const char* lbl, float& value,
    float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_gldragger(const opengl_window& win, const char* lbl, int& value,
    float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2i& value,
    float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3i& value,
    float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4i& value,
    float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_glcheckbox(const opengl_window& win, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}

bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec3f& value) {
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
bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec4f& value) {
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

bool begin_gltreenode(const opengl_window& win, const char* lbl) {
  return ImGui::TreeNode(lbl);
}

void end_gltreenode(const opengl_window& win) { ImGui::TreePop(); }

bool begin_glselectabletreenode(
    const opengl_window& win, const char* lbl, bool& selected) {
  ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow |
                                  ImGuiTreeNodeFlags_OpenOnDoubleClick;
  if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
  auto open = ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
  if (ImGui::IsItemClicked()) selected = true;
  return open;
}

void begin_glselectabletreeleaf(
    const opengl_window& win, const char* lbl, bool& selected) {
  ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_Leaf |
                                  ImGuiTreeNodeFlags_NoTreePushOnOpen;
  if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
  ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
  if (ImGui::IsItemClicked()) selected = true;
}

bool draw_glcombobox(const opengl_window& win, const char* lbl, int& value,
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

bool draw_glcombobox(const opengl_window& win, const char* lbl, string& value,
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

bool draw_glcombobox(const opengl_window& win, const char* lbl, int& idx,
    int num, const function<const char*(int)>& labels, bool include_null) {
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

void begin_glchild(
    const opengl_window& win, const char* lbl, const vec2i& size) {
  ImGui::PushID(lbl);
  ImGui::BeginChild(lbl, ImVec2(size.x, size.y), false);
}
void end_glchild(const opengl_window& win) {
  ImGui::EndChild();
  ImGui::PopID();
}

void draw_glhistogram(
    const opengl_window& win, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec4f>& values) {
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
void        log_glinfo(const opengl_window& win, const char* msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg, "info");
  _log_mutex.unlock();
}
void log_glinfo(const opengl_window& win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_glerror(const opengl_window& win, const char* msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg, "errn");
  _log_mutex.unlock();
}
void log_glerror(const opengl_window& win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_gllogs(const opengl_window& win) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void draw_gllog(const opengl_window& win) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

}  // namespace yocto
