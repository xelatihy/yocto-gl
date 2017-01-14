//
// YGLFW: helper for GLFW code.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
#include "nuklear/nuklear.h"
#include "nuklear/nuklear_glfw_gl3.h"
#include "nuklear/nuklear_glfw_gl2.h"

#include "../yocto/yocto_math.h"

namespace yui {

//
// initialize glfw
//
GLFWwindow* init_glfw(int width, int height, const std::string& title,
                      bool legacy_gl, void* ctx = nullptr,
                      GLFWcharfun text_callback = nullptr);
//
// Clear glfw
//
void clear_glfw(GLFWwindow* window);

//
// Mouse button
//
int mouse_button(GLFWwindow* window);

//
// Mouse position
//
ym::vec2f mouse_pos(GLFWwindow* window);

//
// Window size
//
ym::vec2i window_size(GLFWwindow* window);

//
// Framebuffer size
//
ym::vec2i framebuffer_size(GLFWwindow* window);

//
// Nuklear
//
nk_context* init_nuklear(GLFWwindow* window, bool legacy_gl);

//
// Nuklear
//
void clear_nuklear(nk_context* ctx, bool legacy_gl);

}  // namespace

#endif
