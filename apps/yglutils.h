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

#ifndef _YGLUTILS_H_
#define _YGLUTILS_H_

#include "../yocto/ygl.h"

#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_ext.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

inline unsigned int make_glprogram(const char* vertex, const char* fragment) {
    unsigned int vid = 0, fid = 0, vao = 0;

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    int errflags;
    char errbuf[10000];

    // create vertex
    vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vid, 1, &vertex, NULL);
    glCompileShader(vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }

    // create fragment
    fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fid, 1, &fragment, NULL);
    glCompileShader(fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }

    // create program
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

    return pid;
}

inline void draw_glimage(ygl::vec2i imsize, const void* data, bool as_float,
    ygl::vec2i winsize, ygl::vec2f imcenter, float imscale) {
    static unsigned int gl_txt = 0, gl_prog = 0, gl_texcoord = 0,
                        gl_triangles = 0;

    // initialization
    if (!gl_txt) {
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
        assert(glGetError() == GL_NO_ERROR);
        glGenTextures(1, &gl_txt);
        glBindTexture(GL_TEXTURE_2D, gl_txt);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imsize.x, imsize.y, 0, GL_RGBA,
            (as_float) ? GL_FLOAT : GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        assert(glGetError() == GL_NO_ERROR);
        glGenBuffers(1, &gl_texcoord);
        glBindBuffer(GL_ARRAY_BUFFER, gl_texcoord);
        auto texcoord = std::vector<ygl::vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
        glBufferData(GL_ARRAY_BUFFER, 4 * 2 * sizeof(float), texcoord.data(),
            GL_STATIC_DRAW);
        assert(glGetError() == GL_NO_ERROR);
        glGenBuffers(1, &gl_triangles);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_triangles);
        auto triangles = std::vector<ygl::vec3i>{{0, 1, 2}, {0, 2, 3}};
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * 2 * sizeof(unsigned int),
            triangles.data(), GL_STATIC_DRAW);
        assert(glGetError() == GL_NO_ERROR);
    }

    // update texture
    assert(glGetError() == GL_NO_ERROR);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl_txt);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, imsize.x, imsize.y, GL_RGBA,
        (as_float) ? GL_FLOAT : GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == GL_NO_ERROR);

    // draw
    glUseProgram(gl_prog);
    assert(glGetError() == GL_NO_ERROR);
    glUniform1i(glGetUniformLocation(gl_prog, "txt"), 0);
    glUniform2f(glGetUniformLocation(gl_prog, "winsize"), winsize.x, winsize.y);
    glUniform2f(glGetUniformLocation(gl_prog, "imsize"), imsize.x, imsize.y);
    glUniform2f(glGetUniformLocation(gl_prog, "imcenter"), imcenter.x, imcenter.y);
    glUniform1f(glGetUniformLocation(gl_prog, "imscale"), imscale);
    assert(glGetError() == GL_NO_ERROR);
    glBindBuffer(GL_ARRAY_BUFFER, gl_texcoord);
    glEnableVertexAttribArray(glGetAttribLocation(gl_prog, "texcoord"));
    glVertexAttribPointer(glGetAttribLocation(gl_prog, "texcoord"), 2, GL_FLOAT,
        false, 0, nullptr);
    assert(glGetError() == GL_NO_ERROR);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_triangles);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
    assert(glGetError() == GL_NO_ERROR);
    glUseProgram(0);
    assert(glGetError() == GL_NO_ERROR);
}

inline GLFWwindow* init_glfw_window(int width, int height, const char* title,
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
    if(!gladLoadGL()) throw std::runtime_error("cannot initialize glad");

    return win;
}

inline void init_imgui_widgets(GLFWwindow* win) {
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

#endif
