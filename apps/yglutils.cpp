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

#include "yglutils.h"
#include <cstdarg>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"

#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"

namespace ygl {

void clear_glframebuffer(const vec4f& color, bool clear_depth) {
    glClearColor(color.x, color.y, color.z, color.w);
    if (clear_depth) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
    } else {
        glClear(GL_COLOR_BUFFER_BIT);
    }
}

void set_glviewport(int x, int y, int w, int h) { glViewport(x, y, w, h); }

void set_glviewport(const vec2i& size) { glViewport(0, 0, size.x, size.y); }

void set_glwireframe(bool enabled) {
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

glprogram make_glprogram(const char* vertex, const char* fragment) {
    auto prog = glprogram();

    assert(glGetError() == GL_NO_ERROR);
    glGenVertexArrays(1, &prog.vao);
    glBindVertexArray(prog.vao);
    assert(glGetError() == GL_NO_ERROR);

    int  errflags;
    char errbuf[10000];

    // create vertex
    assert(glGetError() == GL_NO_ERROR);
    prog.vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(prog.vid, 1, &vertex, NULL);
    glCompileShader(prog.vid);
    glGetShaderiv(prog.vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.vid, 10000, 0, errbuf);
        log_error("shader not compiled with error\n{}", errbuf);
        return {};
    }
    assert(glGetError() == GL_NO_ERROR);

    // create fragment
    assert(glGetError() == GL_NO_ERROR);
    prog.fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(prog.fid, 1, &fragment, NULL);
    glCompileShader(prog.fid);
    glGetShaderiv(prog.fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.fid, 10000, 0, errbuf);
        log_error("shader not compiled with error\n{}", errbuf);
        return {};
    }
    assert(glGetError() == GL_NO_ERROR);

    // create program
    assert(glGetError() == GL_NO_ERROR);
    prog.pid = glCreateProgram();
    glAttachShader(prog.pid, prog.vid);
    glAttachShader(prog.pid, prog.fid);
    glLinkProgram(prog.pid);
    glValidateProgram(prog.pid);
    glGetProgramiv(prog.pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        log_error("program not linked with error\n{}", errbuf);
        return {};
    }
    glGetProgramiv(prog.pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        log_error("program not linked with error\n{}", errbuf);
        return {};
    }
    assert(glGetError() == GL_NO_ERROR);

    return prog;
}

gltexture make_gltexture(
    const image<vec4f>& img, bool as_float, bool linear, bool mipmap) {
    auto texture = gltexture();
    assert(glGetError() == GL_NO_ERROR);
    glGenTextures(1, &texture.tid);
    texture.width  = img.width;
    texture.height = img.height;
    glBindTexture(GL_TEXTURE_2D, texture.tid);
    if (as_float) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, img.width, img.height, 0,
            GL_RGBA, GL_FLOAT, img.pixels.data());
    } else {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0,
            GL_RGBA, GL_FLOAT, img.pixels.data());
    }
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pixels.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
    return texture;
}

void update_gltexture(gltexture& texture, const image<vec4f>& img,
    bool as_float, bool linear, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.tid);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_FLOAT, img.pixels.data());
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pixels.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
}

gltexture make_gltexture(
    const image<vec4b>& img, bool as_srgb, bool linear, bool mipmap) {
    auto texture = gltexture();
    assert(glGetError() == GL_NO_ERROR);
    glGenTextures(1, &texture.tid);
    glBindTexture(GL_TEXTURE_2D, texture.tid);
    texture.width  = img.width;
    texture.height = img.height;
    if (as_srgb) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, img.width, img.height, 0,
            GL_RGBA, GL_UNSIGNED_BYTE, img.pixels.data());
    } else {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0,
            GL_RGBA, GL_UNSIGNED_BYTE, img.pixels.data());
    }
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pixels.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
    return texture;
}

void update_gltexture(gltexture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, texture.tid);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_UNSIGNED_BYTE, img.pixels.data());
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pixels.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
}

template <typename T>
glarraybuffer make_glarraybuffer_impl(const vector<T>& data, bool dynamic) {
    auto buf      = glarraybuffer{};
    buf.num       = data.size();
    buf.elem_size = sizeof(T);
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &buf.bid);
    glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(T), data.data(),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return buf;
}

glarraybuffer make_glarraybuffer(const vector<float>& data, bool dynamic) {
    return make_glarraybuffer_impl(data, dynamic);
}
glarraybuffer make_glarraybuffer(const vector<vec2f>& data, bool dynamic) {
    return make_glarraybuffer_impl(data, dynamic);
}
glarraybuffer make_glarraybuffer(const vector<vec3f>& data, bool dynamic) {
    return make_glarraybuffer_impl(data, dynamic);
}
glarraybuffer make_glarraybuffer(const vector<vec4f>& data, bool dynamic) {
    return make_glarraybuffer_impl(data, dynamic);
}

template <typename T>
glelementbuffer make_glelementbuffer_impl(const vector<T>& data, bool dynamic) {
    auto buf      = glelementbuffer{};
    buf.num       = data.size();
    buf.elem_size = sizeof(T);
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &buf.bid);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, data.size() * sizeof(T), data.data(),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return buf;
}

glelementbuffer make_glelementbuffer(const vector<int>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}
glelementbuffer make_glelementbuffer(const vector<vec2i>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}
glelementbuffer make_glelementbuffer(const vector<vec3i>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}

void bind_glprogram(glprogram& prog) { glUseProgram(prog.pid); }
void unbind_glprogram() { glUseProgram(0); }

int get_gluniform_location(const glprogram& prog, const char* name) {
    return glGetUniformLocation(prog.pid, name);
}

void set_gluniform(int loc, int val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform1i(loc, val);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec2i& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform2i(loc, val.x, val.y);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec3i& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform3i(loc, val.x, val.y, val.z);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec4i& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform4i(loc, val.x, val.y, val.z, val.w);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, float val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform1f(loc, val);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec2f& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform2f(loc, val.x, val.y);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec3f& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform3f(loc, val.x, val.y, val.z);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const vec4f& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniform4f(loc, val.x, val.y, val.z, val.w);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const mat4f& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniformMatrix4fv(loc, 1, false, &val.x.x);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int loc, const frame3f& val) {
    assert(glGetError() == GL_NO_ERROR);
    glUniformMatrix4x3fv(loc, 1, false, &val.x.x);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(int loc, const gltexture& texture, int unit) {
    assert(glGetError() == GL_NO_ERROR);
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.tid);
    glUniform1i(loc, unit);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(
    glprogram& prog, const char* var, const gltexture& texture, int unit) {
    set_gluniform_texture(get_gluniform_location(prog, var), texture, unit);
}

void set_gluniform_texture(
    int loc, int loc_on, const gltexture& texture, int unit) {
    assert(glGetError() == GL_NO_ERROR);
    if (texture.tid) {
        glActiveTexture(GL_TEXTURE0 + unit);
        glBindTexture(GL_TEXTURE_2D, texture.tid);
        glUniform1i(loc, unit);
        glUniform1i(loc_on, 1);
    } else {
        glUniform1i(loc_on, 0);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(glprogram& prog, const char* var, const char* var_on,
    const gltexture& texture, int unit) {
    set_gluniform_texture(get_gluniform_location(prog, var),
        get_gluniform_location(prog, var_on), texture, unit);
}

int get_glvertexattrib_location(const glprogram& prog, const char* name) {
    return glGetAttribLocation(prog.pid, name);
}

void set_glvertexattrib(int loc, const glarraybuffer& buf, float val) {
    assert(glGetError() == GL_NO_ERROR);
    if (buf.bid) {
        glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 1, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib1f(loc, val);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec2f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (buf.bid) {
        glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 2, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib2f(loc, val.x, val.y);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec3f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (buf.bid) {
        glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 3, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib3f(loc, val.x, val.y, val.z);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec4f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (buf.bid) {
        glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 4, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib4f(loc, val.x, val.y, val.z, val.w);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void draw_glpoints(const glelementbuffer& buf, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glDrawElements(GL_POINTS, num, GL_UNSIGNED_INT, nullptr);
}

void draw_gllines(const glelementbuffer& buf, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glDrawElements(GL_LINES, num * 2, GL_UNSIGNED_INT, nullptr);
}

void draw_gltriangles(const glelementbuffer& buf, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

void draw_glimage(const gltexture& gl_txt, vec2i imsize, vec2i winsize,
    vec2f imcenter, float imscale) {
    static glprogram       gl_prog      = {};
    static glarraybuffer   gl_texcoord  = {};
    static glelementbuffer gl_triangles = {};

    // initialization
    if (!gl_prog) {
        auto vert   = R"(
            #version 330
            in vec2 texcoord;
            out vec2 frag_texcoord;
            uniform vec2 winsize, imsize;
            uniform vec2 imcenter;
            uniform float imscale;
            void main() {
                vec2 pos = (texcoord - vec2(0.5,0.5)) * imsize * imscale + imcenter;
                gl_Position = vec4(2 * pos.x / winsize.x - 1, 1 - 2 * pos.y / winsize.y, 0, 1);
                frag_texcoord = texcoord;
            }
        )";
        auto frag   = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform sampler2D txt;
            void main() {
                frag_color = texture(txt, frag_texcoord);
            }
        )";
        gl_prog     = make_glprogram(vert, frag);
        gl_texcoord = make_glarraybuffer(
            vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
        gl_triangles = make_glelementbuffer(
            vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    }

    // draw
    bind_glprogram(gl_prog);
    set_gluniform_texture(gl_prog, "txt", gl_txt, 0);
    set_gluniform(gl_prog, "winsize", vec2f{(float)winsize.x, (float)winsize.y});
    set_gluniform(gl_prog, "imsize", vec2f{(float)imsize.x, (float)imsize.y});
    set_gluniform(gl_prog, "imcenter", imcenter);
    set_gluniform(gl_prog, "imscale", imscale);
    set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
    draw_gltriangles(gl_triangles, 2);
    unbind_glprogram();
}

void _glfw_refresh_callback(GLFWwindow* glfw) {
    auto& win = *(const glwindow*)glfwGetWindowUserPointer(glfw);
    if (win.refresh_cb) win.refresh_cb(win);
}

void _glfw_drop_callback(GLFWwindow* glfw, int num, const char** paths) {
    auto& win = *(const glwindow*)glfwGetWindowUserPointer(glfw);
    if (win.drop_cb) {
        auto pathv = vector<string>();
        for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
        win.drop_cb(win, pathv);
    }
}

bool init_glwindow(glwindow& win, int width, int height, const char* title,
    void* user_pointer, std::function<void(const glwindow&)> refresh_cb) {
    // init glfw
    if (!glfwInit()) log_fatal("cannot initialize windowing system");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // create window
    win     = glwindow();
    win.win = glfwCreateWindow(width, height, title, nullptr, nullptr);
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

void delete_glwindow(glwindow& win) {
    glfwDestroyWindow(win.win);
    glfwTerminate();
    win.win = nullptr;
}

void* get_user_pointer(const glwindow& win) { return win.user_ptr; }

void set_drop_callback(glwindow&                                     win,
    function<void(const glwindow& win, const vector<string>& paths)> drop_cb) {
    win.drop_cb = drop_cb;
    glfwSetDropCallback(win.win, _glfw_drop_callback);
}

vec2i get_glframebuffer_size(const glwindow& win) {
    auto size = zero2i;
    glfwGetFramebufferSize(win.win, &size.x, &size.y);
    return size;
}

vec2i get_glwindow_size(const glwindow& win) {
    auto size = zero2i;
    glfwGetWindowSize(win.win, &size.x, &size.y);
    return size;
}

bool should_glwindow_close(const glwindow& win) {
    return glfwWindowShouldClose(win.win);
}

vec2f get_glmouse_pos(const glwindow& win) {
    double mouse_posx, mouse_posy;
    glfwGetCursorPos(win.win, &mouse_posx, &mouse_posy);
    return vec2f{(float)mouse_posx, (float)mouse_posy};
}

bool get_glmouse_left(const glwindow& win) {
    return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
}
bool get_glmouse_right(const glwindow& win) {
    return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
}

bool get_glalt_key(const glwindow& win) {
    return glfwGetKey(win.win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win.win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

bool get_glshift_key(const glwindow& win) {
    return glfwGetKey(win.win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win.win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

void process_glevents(const glwindow& win, bool wait) {
    if (wait)
        glfwWaitEvents();
    else
        glfwPollEvents();
}

void swap_glbuffers(const glwindow& win) { glfwSwapBuffers(win.win); }

void init_glwidgets(const glwindow& win) {
    // init widgets
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename = nullptr;
    ImGui_ImplGlfw_InitForOpenGL(win.win, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    ImGui::StyleColorsDark();
}

bool get_glwidgets_active(const glwindow& win) {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

void begin_glwidgets_frame(const glwindow& win) {
    static auto first_time = true;
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
}

void end_glwidgets_frame(const glwindow& win) {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool begin_glwidgets_window(const glwindow& win, const char* title) {
    return ImGui::Begin(title);
}

bool begin_header_glwidget(const glwindow& win, const char* lbl) {
    if (!ImGui::CollapsingHeader(lbl)) return false;
    ImGui::PushID(lbl);
    return true;
}
void end_header_glwidget(const glwindow& win) { ImGui::PopID(); }

bool draw_button_glwidget(const glwindow& win, const char* lbl) {
    return ImGui::Button(lbl);
}

void draw_label_glwidgets(
    const glwindow& win, const char* lbl, const string& texture) {
    ImGui::LabelText(lbl, "%s", texture.c_str());
}

void draw_label_glwidgets(
    const glwindow& win, const char* lbl, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ImGui::LabelTextV(lbl, fmt, args);
    va_end(args);
}

void draw_separator_glwidget(const glwindow& win) { ImGui::Separator(); }

void continue_glwidgets_line(const glwindow& win) { ImGui::SameLine(); }

bool draw_textinput_glwidget(const glwindow& win, const char* lbl, string& val) {
    char buf[4096];
    auto num = 0;
    for (auto c : val) buf[num++] = c;
    buf[num]    = 0;
    auto edited = ImGui::InputText(lbl, buf, sizeof(buf));
    if (edited) val = buf;
    return edited;
}

bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, float& val, float min, float max) {
    return ImGui::SliderFloat(lbl, &val, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec2f& val, float min, float max) {
    return ImGui::SliderFloat2(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec3f& val, float min, float max) {
    return ImGui::SliderFloat3(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec4f& val, float min, float max) {
    return ImGui::SliderFloat4(lbl, &val.x, min, max);
}

bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, int& val, int min, int max) {
    return ImGui::SliderInt(lbl, &val, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec2i& val, int min, int max) {
    return ImGui::SliderInt2(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec3i& val, int min, int max) {
    return ImGui::SliderInt3(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec4i& val, int min, int max) {
    return ImGui::SliderInt4(lbl, &val.x, min, max);
}

bool draw_dragger_glwidget(const glwindow& win, const char* lbl, float& val,
    float speed, float min, float max) {
    return ImGui::DragFloat(lbl, &val, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec2f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat2(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec3f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat3(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec4f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat4(lbl, &val.x, speed, min, max);
}

bool draw_dragger_glwidget(const glwindow& win, const char* lbl, int& val,
    float speed, int min, int max) {
    return ImGui::DragInt(lbl, &val, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec2i& val,
    float speed, int min, int max) {
    return ImGui::DragInt2(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec3i& val,
    float speed, int min, int max) {
    return ImGui::DragInt3(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec4i& val,
    float speed, int min, int max) {
    return ImGui::DragInt4(lbl, &val.x, speed, min, max);
}

bool draw_checkbox_glwidget(const glwindow& win, const char* lbl, bool& val) {
    return ImGui::Checkbox(lbl, &val);
}

bool draw_coloredit_glwidget(const glwindow& win, const char* lbl, vec3f& val) {
    return ImGui::ColorEdit3(lbl, &val.x);
}

bool draw_coloredit_glwidget(const glwindow& win, const char* lbl, vec4f& val) {
    return ImGui::ColorEdit4(lbl, &val.x);
}

bool begin_treenode_glwidget(const glwindow& win, const char* lbl) {
    return ImGui::TreeNode(lbl);
}

void end_treenode_glwidget(const glwindow& win) { ImGui::TreePop(); }

bool begin_selectabletreenode_glwidget(
    const glwindow& win, const char* lbl, bool& selected) {
    ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow |
                                    ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selected = true;
    return open;
}

void begin_selectabletreeleaf_glwidget(
    const glwindow& win, const char* lbl, bool& selected) {
    ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_Leaf |
                                    ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selected) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(lbl, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selected = true;
}

bool draw_combobox_glwidget(const glwindow& win, const char* lbl, int& val,
    const vector<string>& labels) {
    if (!ImGui::BeginCombo(lbl, labels[val].c_str())) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(labels[i].c_str(), val == i)) val = i;
        if (val == i) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

bool draw_combobox_glwidget(const glwindow& win, const char* lbl, string& val,
    const vector<string>& labels) {
    if (!ImGui::BeginCombo(lbl, val.c_str())) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(labels[i].c_str(), val == labels[i]))
            val = labels[i];
        if (val == labels[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

bool draw_combobox_glwidget(const glwindow& win, const char* lbl, int& idx,
    int num, const function<const char*(int)>& labels, bool include_null) {
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

void begin_child_glwidget(
    const glwindow& win, const char* lbl, const vec2i& size) {
    ImGui::PushID(lbl);
    ImGui::BeginChild(lbl, ImVec2(size.x, size.y), false);
}
void end_child_glwidget(const glwindow& win) {
    ImGui::EndChild();
    ImGui::PopID();
}

}  // namespace ygl
