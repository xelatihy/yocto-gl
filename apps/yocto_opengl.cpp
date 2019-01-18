//
// Utilities to use OpenGL 3, GLFW and ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
#include <cstdarg>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"

#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"

namespace yocto {

void check_opengl_error() {
    if (glGetError() != GL_NO_ERROR) log_error("gl error");
}

void clear_opengl_lframebuffer(const vec4f& color, bool clear_depth) {
    glClearColor(color.x, color.y, color.z, color.w);
    if (clear_depth) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
    } else {
        glClear(GL_COLOR_BUFFER_BIT);
    }
}

void set_opengl_viewport(int x, int y, int w, int h) { glViewport(x, y, w, h); }

void set_opengl_viewport(const vec2i& size) {
    glViewport(0, 0, size.x, size.y);
}

void set_opengl_wireframe(bool enabled) {
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void set_opengl_blending(bool enabled) {
    if (enabled) {
        glEnable(GL_BLEND);
        glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
        glBlendFuncSeparate(
            GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
    } else {
        glDisable(GL_BLEND);
    }
}

bool init_opengl_program(
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
        log_error("shader not compiled with error\n{}", errbuf);
        return false;
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
        log_error("shader not compiled with error\n{}", errbuf);
        return false;
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
        log_error("program not linked with error\n{}", errbuf);
        return false;
    }
    glGetProgramiv(program.program_id, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
        log_error("program not linked with error\n{}", errbuf);
        return false;
    }
    assert(glGetError() == GL_NO_ERROR);
    return true;
}

void delete_opengl_program(opengl_program& program) {
    if (!program) return;
    glDeleteProgram(program.program_id);
    glDeleteShader(program.vertex_shader_id);
    glDeleteShader(program.fragment_shader_id);
    program = {};
}

bool init_opengl_texture(opengl_texture& texture, int width, int height,
    bool as_float, bool as_srgb, bool linear, bool mipmap) {
    texture = opengl_texture();
    assert(glGetError() == GL_NO_ERROR);
    glGenTextures(1, &texture.texture_id);
    texture.width  = width;
    texture.height = height;
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    if (as_float) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA,
            GL_FLOAT, nullptr);
    } else if (as_srgb) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, width, height, 0, GL_RGBA,
            GL_UNSIGNED_BYTE, nullptr);
    } else {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
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
    return true;
}

void update_opengl_texture(
    opengl_texture& texture, const image4f& img, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_FLOAT, data(img));
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == GL_NO_ERROR);
}

void update_opengl_texture_region(opengl_texture& texture, const image4f& img,
    const image_region& region, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    auto clipped = get_image_region(img, region);
    glTexSubImage2D(GL_TEXTURE_2D, 0, region.offsetx, region.offsety,
        region.width, region.height, GL_RGBA, GL_FLOAT, data(clipped));
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == GL_NO_ERROR);
}

void update_opengl_texture(
    opengl_texture& texture, const image4b& img, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_UNSIGNED_BYTE, data(img));
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == GL_NO_ERROR);
}

void update_opengl_texture_region(opengl_texture& texture, const image4b& img,
    const image_region& region, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    auto clipped = get_image_region(img, region);
    glTexSubImage2D(GL_TEXTURE_2D, 0, region.offsetx, region.offsety,
        region.width, region.height, GL_RGBA, GL_UNSIGNED_BYTE, data(clipped));
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == GL_NO_ERROR);
}

void delete_opengl_texture(opengl_texture& texture) {
    if (!texture) return;
    glDeleteTextures(1, &texture.texture_id);
    texture = {};
}

template <typename T>
bool init_opengl_array_buffer_impl(
    opengl_array_buffer& buffer, const vector<T>& array, bool dynamic) {
    buffer           = opengl_array_buffer{};
    buffer.num       = size(array);
    buffer.elem_size = sizeof(T);
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &buffer.buffer_id);
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glBufferData(GL_ARRAY_BUFFER, size(array) * sizeof(T), data(array),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return true;
}

bool init_opengl_array_buffer(
    opengl_array_buffer& buffer, const vector<float>& data, bool dynamic) {
    return init_opengl_array_buffer_impl(buffer, data, dynamic);
}
bool init_opengl_array_buffer(
    opengl_array_buffer& buffer, const vector<vec2f>& data, bool dynamic) {
    return init_opengl_array_buffer_impl(buffer, data, dynamic);
}
bool init_opengl_array_buffer(
    opengl_array_buffer& buffer, const vector<vec3f>& data, bool dynamic) {
    return init_opengl_array_buffer_impl(buffer, data, dynamic);
}
bool init_opengl_array_buffer(
    opengl_array_buffer& buffer, const vector<vec4f>& data, bool dynamic) {
    return init_opengl_array_buffer_impl(buffer, data, dynamic);
}

void delete_opengl_array_buffer(opengl_array_buffer& buffer) {
    if (!buffer) return;
    glDeleteBuffers(1, &buffer.buffer_id);
    buffer = {};
}

template <typename T>
bool init_opengl_elementbuffer_impl(
    opengl_elementbuffer& buffer, const vector<T>& array, bool dynamic) {
    buffer           = opengl_elementbuffer{};
    buffer.num       = size(array);
    buffer.elem_size = sizeof(T);
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &buffer.buffer_id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size(array) * sizeof(T), data(array),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return true;
}

bool init_opengl_elementbuffer(
    opengl_elementbuffer& buffer, const vector<int>& data, bool dynamic) {
    return init_opengl_elementbuffer_impl(buffer, data, dynamic);
}
bool init_opengl_elementbuffer(
    opengl_elementbuffer& buffer, const vector<vec2i>& data, bool dynamic) {
    return init_opengl_elementbuffer_impl(buffer, data, dynamic);
}
bool init_opengl_elementbuffer(
    opengl_elementbuffer& buffer, const vector<vec3i>& data, bool dynamic) {
    return init_opengl_elementbuffer_impl(buffer, data, dynamic);
}

void delete_opengl_elementbuffer(opengl_elementbuffer& buffer) {
    if (!buffer) return;
    glDeleteBuffers(1, &buffer.buffer_id);
    buffer = {};
}

void bind_opengl_program(opengl_program& program) {
    glUseProgram(program.program_id);
}
void unbind_opengl_program() { glUseProgram(0); }

int get_opengl_uniform_location(
    const opengl_program& program, const char* name) {
    return glGetUniformLocation(program.program_id, name);
}

void set_opengl_uniform(int locatiom, int value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform1i(locatiom, value);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec2i& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform2i(locatiom, value.x, value.y);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec3i& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform3i(locatiom, value.x, value.y, value.z);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec4i& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform4i(locatiom, value.x, value.y, value.z, value.w);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, float value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform1f(locatiom, value);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec2f& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform2f(locatiom, value.x, value.y);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec3f& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform3f(locatiom, value.x, value.y, value.z);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const vec4f& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform4f(locatiom, value.x, value.y, value.z, value.w);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const mat4f& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniformMatrix4fv(locatiom, 1, false, &value.x.x);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform(int locatiom, const frame3f& value) {
    assert(glGetError() == GL_NO_ERROR);
    glUniformMatrix4x3fv(locatiom, 1, false, &value.x.x);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform_texture(
    int locatiom, const opengl_texture& texture, int unit) {
    assert(glGetError() == GL_NO_ERROR);
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glUniform1i(locatiom, unit);
    assert(glGetError() == GL_NO_ERROR);
}

void set_opengl_uniform_texture(opengl_program& program, const char* name,
    const opengl_texture& texture, int unit) {
    set_opengl_uniform_texture(
        get_opengl_uniform_location(program, name), texture, unit);
}

void set_opengl_uniform_texture(
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

void set_opengl_uniform_texture(opengl_program& program, const char* name,
    const char* name_on, const opengl_texture& texture, int unit) {
    set_opengl_uniform_texture(get_opengl_uniform_location(program, name),
        get_opengl_uniform_location(program, name_on), texture, unit);
}

int get_opengl_vertexattrib_location(
    const opengl_program& program, const char* name) {
    return glGetAttribLocation(program.program_id, name);
}

void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, float value) {
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

void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec2f& value) {
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

void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec3f& value) {
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

void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec4f& value) {
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

void draw_opengl_points(const opengl_elementbuffer& buffer, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
    glDrawElements(GL_POINTS, num, GL_UNSIGNED_INT, nullptr);
}

void draw_opengl_lines(const opengl_elementbuffer& buffer, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
    glDrawElements(GL_LINES, num * 2, GL_UNSIGNED_INT, nullptr);
}

void draw_opengl_triangles(const opengl_elementbuffer& buffer, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
    glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

void draw_opengl_image(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale) {
    static opengl_program       gl_prog      = {};
    static opengl_array_buffer  gl_texcoord  = {};
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
        init_opengl_program(gl_prog, vert, frag);
        init_opengl_array_buffer(
            gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
        init_opengl_elementbuffer(
            gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    }

    // draw
    check_opengl_error();
    bind_opengl_program(gl_prog);
    set_opengl_uniform_texture(gl_prog, "txt", texture, 0);
    set_opengl_uniform(
        gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
    set_opengl_uniform(gl_prog, "image_size",
        vec2f{(float)texture.width, (float)texture.height});
    set_opengl_uniform(gl_prog, "image_center", image_center);
    set_opengl_uniform(gl_prog, "image_scale", image_scale);
    set_opengl_vertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
    draw_opengl_triangles(gl_triangles, 2);
    unbind_opengl_program();
    check_opengl_error();
}

void draw_opengl_image_background(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size) {
    static opengl_program       gl_prog      = {};
    static opengl_array_buffer  gl_texcoord  = {};
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
        init_opengl_program(gl_prog, vert, frag);
        init_opengl_array_buffer(
            gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
        init_opengl_elementbuffer(
            gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    }

    // draw
    bind_opengl_program(gl_prog);
    set_opengl_uniform(
        gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
    set_opengl_uniform(gl_prog, "image_size",
        vec2f{(float)texture.width, (float)texture.height});
    set_opengl_uniform(
        gl_prog, "border_size", vec2f{(float)border_size, (float)border_size});
    set_opengl_uniform(gl_prog, "image_center", image_center);
    set_opengl_uniform(gl_prog, "image_scale", image_scale);
    set_opengl_vertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
    draw_opengl_triangles(gl_triangles, 2);
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

bool init_opengl_window(opengl_window& win, const vec2i& size,
    const string& title, void* user_pointer,
    refresh_opengl_callback refresh_cb) {
    // init glfw
    if (!glfwInit()) log_fatal("cannot initialize windowing system");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // create window
    win     = opengl_window();
    win.win = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
    if (!win.win) return false;
    glfwMakeContextCurrent(win.win);
    glfwSwapInterval(1);  // Enable vsync

    // set user data
    glfwSetWindowRefreshCallback(win.win, _glfw_refresh_callback);
    glfwSetWindowUserPointer(win.win, &win);
    win.user_ptr   = user_pointer;
    win.refresh_cb = refresh_cb;

    // init gl extensions
    if (!gladLoadGL()) log_fatal("cannot initialize OpenGL extensions");

    return true;
}

void delete_opengl_window(opengl_window& win) {
    glfwDestroyWindow(win.win);
    glfwTerminate();
    win.win = nullptr;
}

void* get_opengl_user_pointer(const opengl_window& win) { return win.user_ptr; }

void set_drop_opengl_callback(
    opengl_window& win, drop_opengl_callback drop_cb) {
    win.drop_cb = drop_cb;
    glfwSetDropCallback(win.win, _glfw_drop_callback);
}

vec2i get_opengl_framebuffer_size(const opengl_window& win) {
    auto size = zero2i;
    glfwGetFramebufferSize(win.win, &size.x, &size.y);
    return size;
}

vec2i get_opengl_window_size(const opengl_window& win) {
    auto size = zero2i;
    glfwGetWindowSize(win.win, &size.x, &size.y);
    return size;
}

bool should_opengl_window_close(const opengl_window& win) {
    return glfwWindowShouldClose(win.win);
}

vec2f get_opengl_mouse_pos(const opengl_window& win) {
    double mouse_posx, mouse_posy;
    glfwGetCursorPos(win.win, &mouse_posx, &mouse_posy);
    return vec2f{(float)mouse_posx, (float)mouse_posy};
}

bool get_opengl_mouse_left(const opengl_window& win) {
    return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
}
bool get_opengl_mouse_right(const opengl_window& win) {
    return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
}

bool get_opengl_alt_key(const opengl_window& win) {
    return glfwGetKey(win.win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win.win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

bool get_opengl_shift_key(const opengl_window& win) {
    return glfwGetKey(win.win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win.win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

void process_opengl_events(const opengl_window& win, bool wait) {
    if (wait)
        glfwWaitEvents();
    else
        glfwPollEvents();
}

void swap_opengl_buffers(const opengl_window& win) { glfwSwapBuffers(win.win); }

void init_opengl_widgets(const opengl_window& win) {
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
}

bool get_opengl_widgets_active(const opengl_window& win) {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

void begin_opengl_widgets_frame(const opengl_window& win) {
    static auto first_time = true;
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 360});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
}

void end_opengl_widgets_frame(const opengl_window& win) {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool begin_opengl_widgets_window(const opengl_window& win, const char* title) {
    return ImGui::Begin(title);
}

bool begin_header_opengl_widget(const opengl_window& win, const char* lbl) {
    if (!ImGui::CollapsingHeader(lbl)) return false;
    ImGui::PushID(lbl);
    return true;
}
void end_header_opengl_widget(const opengl_window& win) { ImGui::PopID(); }

bool begin_tabbar_opengl_widget(const opengl_window& win, const char* lbl) {
    if (!ImGui::BeginTabBar(lbl)) return false;
    ImGui::PushID(lbl);
    return true;
}
void end_tabbar_opengl_widget(const opengl_window& win) {
    ImGui::PopID();
    ImGui::EndTabBar();
}

bool begin_tabitem_opengl_widget(const opengl_window& win, const char* lbl) {
    if (!ImGui::BeginTabItem(lbl)) return false;
    ImGui::PushID(lbl);
    return true;
}
void end_tabitem_opengl_widget(const opengl_window& win) {
    ImGui::PopID();
    ImGui::EndTabItem();
}

bool draw_button_opengl_widget(const opengl_window& win, const char* lbl) {
    return ImGui::Button(lbl);
}

void draw_label_opengl_widget(
    const opengl_window& win, const char* lbl, const string& texture) {
    ImGui::LabelText(lbl, "%s", texture.c_str());
}

void draw_label_opengl_widget(
    const opengl_window& win, const char* lbl, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ImGui::LabelTextV(lbl, fmt, args);
    va_end(args);
}

void draw_separator_opengl_widget(const opengl_window& win) {
    ImGui::Separator();
}

void continue_opengl_widget_line(const opengl_window& win) {
    ImGui::SameLine();
}

bool draw_textinput_opengl_widget(
    const opengl_window& win, const char* lbl, string& value) {
    char buffer[4096];
    auto num = 0;
    for (auto c : value) buffer[num++] = c;
    buffer[num] = 0;
    auto edited = ImGui::InputText(lbl, buffer, sizeof(buffer));
    if (edited) value = buffer;
    return edited;
}

bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    float& value, float min, float max) {
    return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec2f& value, float min, float max) {
    return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec3f& value, float min, float max) {
    return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec4f& value, float min, float max) {
    return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, int& value, int min, int max) {
    return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec2i& value, int min, int max) {
    return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec3i& value, int min, int max) {
    return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec4i& value, int min, int max) {
    return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    float& value, float speed, float min, float max) {
    return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec1f& value, float speed, float min, float max) {
    return ImGui::DragFloat(lbl, &value.x, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec2f& value, float speed, float min, float max) {
    return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec3f& value, float speed, float min, float max) {
    return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec4f& value, float speed, float min, float max) {
    return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    int& value, float speed, int min, int max) {
    return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec2i& value, float speed, int min, int max) {
    return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec3i& value, float speed, int min, int max) {
    return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec4i& value, float speed, int min, int max) {
    return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_checkbox_opengl_widget(
    const opengl_window& win, const char* lbl, bool& value) {
    return ImGui::Checkbox(lbl, &value);
}

bool draw_coloredit_opengl_widget(
    const opengl_window& win, const char* lbl, vec3f& value) {
    return ImGui::ColorEdit3(lbl, &value.x);
}

bool draw_coloredit_opengl_widget(
    const opengl_window& win, const char* lbl, vec4f& value) {
    return ImGui::ColorEdit4(lbl, &value.x);
}

bool begin_treenode_opengl_widget(const opengl_window& win, const char* lbl) {
    return ImGui::TreeNode(lbl);
}

void end_treenode_opengl_widget(const opengl_window& win) { ImGui::TreePop(); }

bool begin_selectabletreenode_opengl_widget(
    const opengl_window& win, const char* lbl, bool& selected) {
    ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow |
                                    ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selected = true;
    return open;
}

void begin_selectabletreeleaf_opengl_widget(
    const opengl_window& win, const char* lbl, bool& selected) {
    ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_Leaf |
                                    ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selected = true;
}

bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    int& value, const vector<string>& labels) {
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

bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    string& value, const vector<string>& labels) {
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

bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    int& idx, int num, const function<const char*(int)>& labels,
    bool include_null) {
    if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx) : "<none>"))
        return false;
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

void begin_child_opengl_widget(
    const opengl_window& win, const char* lbl, const vec2i& size) {
    ImGui::PushID(lbl);
    ImGui::BeginChild(lbl, ImVec2(size.x, size.y), false);
}
void end_child_opengl_widget(const opengl_window& win) {
    ImGui::EndChild();
    ImGui::PopID();
}

}  // namespace yocto
