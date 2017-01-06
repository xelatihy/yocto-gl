//
// YGLFW: helper for GLFW code.
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
#include <OpenGL/gl3.h>
#endif
#include <GLFW/glfw3.h>
// clang-format on

#include <cstring>
#include <string>

#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_VERTEX_BUFFER_OUTPUT
#define NK_INCLUDE_FONT_BAKING
#define NK_INCLUDE_DEFAULT_FONT
#define NK_IMPLEMENTATION
#define NK_GLFW_GL3_IMPLEMENTATION
#define NK_GLFW_GL2_IMPLEMENTATION
#include "nuklear/nuklear.h"
#include "nuklear/nuklear_glfw_gl3.h"
#include "nuklear/nuklear_glfw_gl2.h"

#include "../yocto/yocto_math.h"

namespace yui {
//
// Support
//
static inline void _glfw_error_callback(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

//
// initialize glfw
//
inline GLFWwindow* init_glfw(int width, int height, const std::string& title,
                             bool legacy_gl, void* ctx = nullptr,
                             GLFWcharfun text_callback = nullptr) {
    // window
    if (!glfwInit()) return nullptr;

    // profile creation
    if (!legacy_gl) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#if __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    }

    auto window = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, ctx);

    glfwSetErrorCallback(_glfw_error_callback);
    if (text_callback) glfwSetCharCallback(window, text_callback);

    return window;
}

//
// Clear glfw
//
inline void clear_glfw(GLFWwindow* window) {
    glfwDestroyWindow(window);
    glfwTerminate();
}

//
// Mouse button
//
inline int mouse_button(GLFWwindow* window) {
    auto mouse1 = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
    if (mouse1) return 1;
    if (mouse2) return 2;
    if (mouse3) return 3;
#if 0
        if (action == GLFW_RELEASE) {
            vparams.mouse_button = 0;
        } else if (button == GLFW_MOUSE_BUTTON_1 && !mods) {
            vparams.mouse_button = 1;
        } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_CONTROL)) {
            vparams.mouse_button = 2;
        } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_SHIFT)) {
            vparams.mouse_button = 3;
        } else if (button == GLFW_MOUSE_BUTTON_2) {
            vparams.mouse_button = 2;
        } else {
            vparams.mouse_button = 0;
        }
#endif
    return 0;
}

//
// Mouse position
//
inline ym::vec2f mouse_pos(GLFWwindow* window) {
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    return {(float)x, (float)y};
}

//
// Window size
//
inline ym::vec2i window_size(GLFWwindow* window) {
    auto ret = ym::zero2i;
    glfwGetWindowSize(window, &ret[0], &ret[1]);
    return ret;
}

//
// Framebuffer size
//
inline ym::vec2i framebuffer_size(GLFWwindow* window) {
    auto ret = ym::zero2i;
    glfwGetFramebufferSize(window, &ret[0], &ret[1]);
    return ret;
}

//
// Nuklear
//
inline nk_context* init_nuklear(GLFWwindow* window, bool legacy_gl) {
    if (legacy_gl) {
        auto nuklear_ctx = nk_glfw3_gl2_init(window, NK_GLFW3_GL2_DEFAULT);
        nk_font_atlas* atlas;
        nk_glfw3_gl2_font_stash_begin(&atlas);
        nk_glfw3_gl2_font_stash_end();
        return nuklear_ctx;
    } else {
        auto nuklear_ctx = nk_glfw3_gl3_init(window, NK_GLFW3_GL3_DEFAULT);
        nk_font_atlas* atlas;
        nk_glfw3_gl3_font_stash_begin(&atlas);
        nk_glfw3_gl3_font_stash_end();
        return nuklear_ctx;
    }
}

//
// Nuklear
//
inline void clear_nuklear(nk_context* ctx, bool legacy_gl) {
    if (legacy_gl) {
    } else {
        nk_glfw3_gl3_shutdown();
    }
}

}  // namespace

#endif
