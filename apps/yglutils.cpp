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
#include "glad/glad.h"

#include <GLFW/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

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

uint make_glprogram(const char* vertex, const char* fragment) {
    uint vid = 0, fid = 0, vao = 0;

    assert(glGetError() == GL_NO_ERROR);
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    assert(glGetError() == GL_NO_ERROR);

    int errflags;
    char errbuf[10000];

    // create vertex
    assert(glGetError() == GL_NO_ERROR);
    vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vid, 1, &vertex, NULL);
    glCompileShader(vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create fragment
    assert(glGetError() == GL_NO_ERROR);
    fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fid, 1, &fragment, NULL);
    glCompileShader(fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create program
    assert(glGetError() == GL_NO_ERROR);
    auto pid = glCreateProgram();
    glAttachShader(pid, vid);
    glAttachShader(pid, fid);
    glLinkProgram(pid);
    glValidateProgram(pid);
    glGetProgramiv(pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    return pid;
}

uint make_gltexture(const image4f& img, bool linear, bool mipmap) {
    auto tid = (uint)0;
    assert(glGetError() == GL_NO_ERROR);
    glGenTextures(1, &tid);
    glBindTexture(GL_TEXTURE_2D, tid);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0, GL_RGBA,
        GL_FLOAT, img.pxl.data());
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pxl.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
    return tid;
}

void update_gltexture(int tid, const image4f& img, bool linear, bool mipmap) {
    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, tid);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_FLOAT, img.pxl.data());
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (!img.pxl.empty()) glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    assert(glGetError() == GL_NO_ERROR);
}

template <typename T>
uint make_glarraybuffer_impl(const std::vector<T>& buf, bool dynamic) {
    auto bid = (uint)0;
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &bid);
    glBindBuffer(GL_ARRAY_BUFFER, bid);
    glBufferData(GL_ARRAY_BUFFER, buf.size() * sizeof(T), buf.data(),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return bid;
}

uint make_glarraybuffer(const std::vector<float>& buf, bool dynamic) {
    return make_glarraybuffer_impl(buf, dynamic);
}
uint make_glarraybuffer(const std::vector<vec2f>& buf, bool dynamic) {
    return make_glarraybuffer_impl(buf, dynamic);
}
uint make_glarraybuffer(const std::vector<vec3f>& buf, bool dynamic) {
    return make_glarraybuffer_impl(buf, dynamic);
}
uint make_glarraybuffer(const std::vector<vec4f>& buf, bool dynamic) {
    return make_glarraybuffer_impl(buf, dynamic);
}

template <typename T>
uint make_glelementbuffer_impl(const std::vector<T>& buf, bool dynamic) {
    auto bid = (uint)0;
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &bid);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, buf.size() * sizeof(T), buf.data(),
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    assert(glGetError() == GL_NO_ERROR);
    return bid;
}

uint make_glelementbuffer(const std::vector<int>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}
uint make_glelementbuffer(const std::vector<vec2i>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}
uint make_glelementbuffer(const std::vector<vec3i>& buf, bool dynamic) {
    return make_glelementbuffer_impl(buf, dynamic);
}

void bind_glprogram(uint pid) { glUseProgram(pid); }

int get_gluniform_location(uint pid, const char* name) {
    return glGetUniformLocation(pid, name);
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

void set_gluniform_texture(int loc, uint tid, int unit) {
    assert(glGetError() == GL_NO_ERROR);
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, tid);
    glUniform1i(loc, unit);
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(uint pid, const char* var, uint tid, int unit) {
    set_gluniform_texture(glGetUniformLocation(pid, var), tid, unit);
}

void set_gluniform_texture(int loc, int loc_on, uint tid, int unit) {
    assert(glGetError() == GL_NO_ERROR);
    if (tid) {
        glActiveTexture(GL_TEXTURE0 + unit);
        glBindTexture(GL_TEXTURE_2D, tid);
        glUniform1i(loc, unit);
        glUniform1i(loc_on, 1);
    } else {
        glUniform1i(loc_on, 0);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(
    uint pid, const char* var, const char* var_on, uint tid, int unit) {
    set_gluniform_texture(glGetUniformLocation(pid, var),
        glGetUniformLocation(pid, var_on), tid, unit);
}

int get_glvertexattrib_location(uint pid, const char* name) {
    return glGetAttribLocation(pid, name);
}

void set_glvertexattrib(int loc, int bid, float val) {
    assert(glGetError() == GL_NO_ERROR);
    if (bid) {
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 1, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib1f(loc, val);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, int bid, const vec2f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (bid) {
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 2, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib2f(loc, val.x, val.y);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, int bid, const vec3f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (bid) {
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 3, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib3f(loc, val.x, val.y, val.z);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(int loc, int bid, const vec4f& val) {
    assert(glGetError() == GL_NO_ERROR);
    if (bid) {
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glEnableVertexAttribArray(loc);
        glVertexAttribPointer(loc, 4, GL_FLOAT, false, 0, nullptr);
    } else {
        glVertexAttrib4f(loc, val.x, val.y, val.z, val.w);
    }
    assert(glGetError() == GL_NO_ERROR);
}

void draw_glpoints(uint bid, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_POINTS, num, GL_UNSIGNED_INT, nullptr);
}

void draw_gllines(uint bid, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_LINES, num * 2, GL_UNSIGNED_INT, nullptr);
}

void draw_gltriangles(uint bid, int num) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

void draw_glimage(
    uint gl_txt, vec2i imsize, vec2i winsize, vec2f imcenter, float imscale) {
    static uint gl_prog = 0, gl_texcoord = 0, gl_triangles = 0;

    // initialization
    if (!gl_prog) {
        auto vert = R"(
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
        auto frag = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform sampler2D txt;
            void main() {
                frag_color = texture(txt, frag_texcoord);
            }
        )";
        gl_prog = make_glprogram(vert, frag);
        gl_texcoord = make_glarraybuffer(
            std::vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
        gl_triangles = make_glelementbuffer(
            std::vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    }

    // draw
    bind_glprogram(gl_prog);
    set_gluniform_texture(gl_prog, "txt", gl_txt, 0);
    set_gluniform(
        gl_prog, "winsize", vec2f{(float)winsize.x, (float)winsize.y});
    set_gluniform(gl_prog, "imsize", vec2f{(float)imsize.x, (float)imsize.y});
    set_gluniform(gl_prog, "imcenter", imcenter);
    set_gluniform(gl_prog, "imscale", imscale);
    set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
    draw_gltriangles(gl_triangles, 2);
    bind_glprogram(0);
}

glwindow* make_glwindow(int width, int height, const char* title,
    void* user_pointer, void (*refresh)(GLFWwindow*)) {
    // init glfw
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // create window
    auto win = glfwCreateWindow(width, height, title, nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    // set user data
    glfwSetWindowRefreshCallback(win, refresh);
    glfwSetWindowUserPointer(win, user_pointer);

    // init gl extensions
    if (!gladLoadGL()) throw std::runtime_error("cannot initialize glad");

    return win;
}

void delete_glwindow(glwindow* win) {
    glfwDestroyWindow(win);
    glfwTerminate();
}

void* get_user_pointer(glwindow* win) { return glfwGetWindowUserPointer(win); }

vec2i get_glframebuffer_size(glwindow* win) {
    auto size = zero2i;
    glfwGetFramebufferSize(win, &size.x, &size.y);
    return size;
}

vec2i get_glwindow_size(glwindow* win) {
    auto size = zero2i;
    glfwGetWindowSize(win, &size.x, &size.y);
    return size;
}

bool should_glwindow_close(glwindow* win) { return glfwWindowShouldClose(win); }

vec2f get_glmouse_pos(glwindow* win) {
    double mouse_posx, mouse_posy;
    glfwGetCursorPos(win, &mouse_posx, &mouse_posy);
    return vec2f{(float)mouse_posx, (float)mouse_posy};
}

bool get_glmouse_left(glwindow* win) {
    return glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
}
bool get_glmouse_right(glwindow* win) {
    return glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
}

bool get_glalt_key(glwindow* win) {
    return glfwGetKey(win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

bool get_glshift_key(glwindow* win) {
    return glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

void process_glevents(glwindow* win, bool wait) {
    if (wait)
        glfwWaitEvents();
    else
        glfwPollEvents();
}

void swap_glbuffers(glwindow* win) { glfwSwapBuffers(win); }

void init_glwidgets(glwindow* win) {
    // init widgets
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename = nullptr;
    ImGui_ImplGlfw_InitForOpenGL(win, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    ImGui::StyleColorsDark();
}

bool get_glwidgets_active(glwindow* win) {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

void begin_glwidgets_frame(glwindow* win) {
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

void end_glwidgets_frame(glwindow* win) {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool begin_glwidgets_window(glwindow* win, const char* title) {
    return ImGui::Begin(title);
}

void draw_imgui_label(glwindow* win, const char* lbl, const std::string& txt) {
    ImGui::LabelText(lbl, "%s", txt.c_str());
}

void draw_imgui_label(glwindow* win, const char* lbl, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ImGui::LabelTextV(lbl, fmt, args);
    va_end(args);
}

void draw_separator_glwidget(glwindow* win) { ImGui::Separator(); }

void continue_glwidgets_line(glwindow* win) { ImGui::SameLine(); }

bool draw_inputtext_glwidget(glwindow* win, const char* lbl, std::string& val) {
    char buf[4096];
    auto num = 0;
    for (auto c : val) buf[num++] = c;
    buf[num] = 0;
    auto edited = ImGui::InputText(lbl, buf, sizeof(buf));
    if (edited) val = buf;
    return edited;
}

bool draw_slider_glwidget(
    glwindow* win, const char* lbl, float& val, float min, float max) {
    return ImGui::SliderFloat(lbl, &val, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec2f& val, float min, float max) {
    return ImGui::SliderFloat2(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec3f& val, float min, float max) {
    return ImGui::SliderFloat3(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec4f& val, float min, float max) {
    return ImGui::SliderFloat4(lbl, &val.x, min, max);
}

bool draw_slider_glwidget(
    glwindow* win, const char* lbl, int& val, int min, int max) {
    return ImGui::SliderInt(lbl, &val, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec2i& val, int min, int max) {
    return ImGui::SliderInt2(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec3i& val, int min, int max) {
    return ImGui::SliderInt3(lbl, &val.x, min, max);
}
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec4i& val, int min, int max) {
    return ImGui::SliderInt4(lbl, &val.x, min, max);
}

bool draw_dragger_glwidget(glwindow* win, const char* lbl, float& val,
    float speed, float min, float max) {
    return ImGui::DragFloat(lbl, &val, speed, min, max);
}
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec2f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat2(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec3f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat3(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec4f& val,
    float speed, float min, float max) {
    return ImGui::DragFloat4(lbl, &val.x, speed, min, max);
}

bool draw_dragger_glwidget(
    glwindow* win, const char* lbl, int& val, float speed, int min, int max) {
    return ImGui::DragInt(lbl, &val, speed, min, max);
}
bool draw_dragger_glwidget(
    glwindow* win, const char* lbl, vec2i& val, float speed, int min, int max) {
    return ImGui::DragInt2(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(
    glwindow* win, const char* lbl, vec3i& val, float speed, int min, int max) {
    return ImGui::DragInt3(lbl, &val.x, speed, min, max);
}
bool draw_dragger_glwidget(
    glwindow* win, const char* lbl, vec4i& val, float speed, int min, int max) {
    return ImGui::DragInt4(lbl, &val.x, speed, min, max);
}

bool draw_checkbox_glwidget(glwindow* win, const char* lbl, bool& val) {
    return ImGui::Checkbox(lbl, &val);
}

bool draw_coloredit_glwidget(glwindow* win, const char* lbl, vec3f& val) {
    return ImGui::ColorEdit3(lbl, &val.x);
}

bool draw_coloredit_glwidget(glwindow* win, const char* lbl, vec4f& val) {
    return ImGui::ColorEdit4(lbl, &val.x);
}

bool begin_treenode_glwidget(glwindow* win, const char* lbl) {
    return ImGui::TreeNode(lbl);
}

void end_treenode_glwidget(glwindow* win) { ImGui::TreePop(); }

bool begin_selectabletreenode_glwidget(
    glwindow* win, const char* lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selection = content;
    return open;
}

void begin_selectabletreeleaf_glwidget(
    glwindow* win, const char* lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selection = content;
}

bool draw_combobox_glwidget(glwindow* win, const char* lbl, int& val,
    const std::vector<std::string>& labels) {
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

bool draw_combobox_glwidget(glwindow* win, const char* lbl, std::string& val,
    const std::vector<std::string>& labels) {
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

bool draw_combobox_glwidget(glwindow* win, const char* lbl, int& idx,
    const std::vector<void*>& vals, const char* (*label)(void*)) {
    if (!ImGui::BeginCombo(lbl, label(vals.at(idx)))) return false;
    auto old_idx = idx;
    for (auto i = 0; i < vals.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(label(vals.at(i)), idx == i)) idx = i;
        if (idx == i) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return idx != old_idx;
}

bool draw_combobox_glwidget(glwindow* win, const char* lbl, void*& val,
    const std::vector<void*>& vals, const char* (*label)(void*),
    bool include_null) {
    if (!ImGui::BeginCombo(lbl, (val) ? label(val) : "<none>")) return false;
    auto old_val = val;
    if (include_null) {
        ImGui::PushID(100000);
        if (ImGui::Selectable("<none>", val == nullptr)) val = nullptr;
        if (val == nullptr) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    for (auto i = 0; i < vals.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(label(vals[i]), val == vals[i])) val = vals[i];
        if (val == vals[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

void begin_child_glwidget(glwindow* win, const char* lbl, const vec2i& size) {
    ImGui::PushID(lbl);
    ImGui::BeginChild(lbl, ImVec2(size.x, size.y), false);
}
void end_child_glwidget(glwindow* win) {
    ImGui::EndChild();
    ImGui::PopID();
}

}  // namespace ygl
