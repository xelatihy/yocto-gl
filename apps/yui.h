//
// YUI: helpers to define a user interface with GLFW.
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _YUI_H_
#define _YUI_H_

// clang-format off
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl.h>
#endif
#include <GLFW/glfw3.h>
// clang-format on

#include "../yocto/yocto_math.h"
#include <functional>
#include <string>

namespace yui {

struct info {
    ym::vec2i win_size;                 // window size
    ym::vec2i framebuffer_size;         // framebuffer size
    int mouse_button = 0;               // mouse button
    ym::vec2f mouse_pos = ym::zero2f;   // mouse position
    ym::vec2f mouse_last = ym::zero2f;  // last mouse position
};

using init_cb = std::function<void(const info& info)>;
using text_cb = std::function<void(const info& info, unsigned int key)>;
using window_size_cb = std::function<void(const info& info)>;
using framebuffer_size_cb = std::function<void(const info& info)>;
using mouse_button_cb = std::function<void(const info& info)>;
using mouse_pos_cb = std::function<void(const info& info)>;
using window_refresh_cb = std::function<void(const info& info)>;
using update_cb = std::function<int(const info& info)>;

struct context {
    // configuring options
    std::vector<init_cb> init;
    std::vector<text_cb> text;
    std::vector<window_size_cb> window_size;
    std::vector<framebuffer_size_cb> framebuffer_size;
    std::vector<mouse_button_cb> mouse_button;
    std::vector<mouse_pos_cb> mouse_pos;
    std::vector<window_refresh_cb> window_refresh;
    std::vector<update_cb> update;

    // information
    info _info;  // current info
};

// callbacks
static inline void _init_callback(GLFWwindow* window) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    for (auto& cb : ctx.init) cb(ctx._info);
}

static inline void _text_callback(GLFWwindow* window, unsigned int key) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    for (auto& cb : ctx.text) cb(ctx._info, key);
}

static inline void _window_size_callback(GLFWwindow* window, int w, int h) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    ctx._info.win_size = {w, h};
    for (auto& cb : ctx.window_size) cb(ctx._info);
}

static inline void _framebuffer_size_callback(GLFWwindow* window, int w,
                                              int h) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    glViewport(0, 0, w, h);
    ctx._info.framebuffer_size = {w, h};
    for (auto& cb : ctx.framebuffer_size) cb(ctx._info);
}

static inline void _mouse_button_callback(GLFWwindow* window, int button,
                                          int action, int mods) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    if (action == GLFW_RELEASE) {
        ctx._info.mouse_button = 0;
    } else if (button == GLFW_MOUSE_BUTTON_1 && !mods) {
        ctx._info.mouse_button = 1;
    } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_CONTROL)) {
        ctx._info.mouse_button = 2;
    } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_SHIFT)) {
        ctx._info.mouse_button = 3;
    } else if (button == GLFW_MOUSE_BUTTON_2) {
        ctx._info.mouse_button = 2;
    } else {
        ctx._info.mouse_button = 0;
    }
    for (auto& cb : ctx.mouse_button) cb(ctx._info);
}

static inline void _mouse_pos_callback(GLFWwindow* window, double x, double y) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    ctx._info.mouse_last = ctx._info.mouse_pos;
    ctx._info.mouse_pos = {(float)x, (float)y};
    for (auto& cb : ctx.mouse_pos) cb(ctx._info);
}

static inline void _window_refresh_callback(GLFWwindow* window) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    for (auto& cb : ctx.window_refresh) cb(ctx._info);
    glfwSwapBuffers(window);
}

static inline int _update_callback(GLFWwindow* window) {
    auto& ctx = *(context*)glfwGetWindowUserPointer(window);
    auto action = 0;
    for (auto& cb : ctx.update) {
        auto up = cb(ctx._info);
        if (up < 0)
            action = -1;
        else
            action = ym::max(action, up);
    }
    return action;
}

static inline void _error_callback(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

static inline void _print_gamma_ramp() {
    auto& ramp = *glfwGetGammaRamp(glfwGetPrimaryMonitor());
    for (auto i = 0; i < ramp.size; i++) {
        printf("%03d %3.3f %3.3f %3.3f\n", i, ramp.red[i] / 65535.0f,
               ramp.green[i] / 65535.0f, ramp.blue[i] / 65535.0f);
    }
}

inline bool ui_loop(context& context, int width, int height,
                    const std::string& title, bool modern = false) {
    // window
    if (!glfwInit()) return false;

    // profile creation
    if (modern) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    }

    auto window = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, &context);

    // callbacks
    glfwSetErrorCallback(_error_callback);
    glfwSetCharCallback(window, _text_callback);
    glfwSetWindowSizeCallback(window, _window_size_callback);
    glfwSetFramebufferSizeCallback(window, _framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, _mouse_button_callback);
    glfwSetCursorPosCallback(window, _mouse_pos_callback);
    glfwSetWindowRefreshCallback(window, _window_refresh_callback);

    // get initial values
    int w, h;
    glfwGetWindowSize(window, &w, &h);
    context._info.win_size = {w, h};
    glfwGetFramebufferSize(window, &w, &h);
    context._info.framebuffer_size = {w, h};

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) exit(EXIT_FAILURE);
#endif

    // load textures
    _init_callback(window);

    // ui loop
    while (!glfwWindowShouldClose(window)) {
        _window_refresh_callback(window);
        auto update = _update_callback(window);
        if (update < 0) break;
        if (update > 0)
            glfwPollEvents();
        else
            glfwWaitEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    return true;
}

}  // namespace

#endif
