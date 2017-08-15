//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#include "yocto_glu.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>

// clang-format off
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl3.h>
#endif
// clang-format on

namespace yglu {

//
// Checks for GL error and then prints
//
bool check_error(bool print) {
    auto ok = glGetError();
    if (ok == GL_NO_ERROR) return true;
    if (!print) return false;
    switch (ok) {
        case GL_NO_ERROR: printf("GL_NO_ERROR\n"); break;
        case GL_INVALID_ENUM: printf("GL_INVALID_ENUM\n"); break;
        case GL_INVALID_VALUE: printf("GL_INVALID_VALUE\n"); break;
        case GL_INVALID_OPERATION: printf("GL_INVALID_OPERATION\n"); break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            printf("GL_INVALID_FRAMEBUFFER_OPERATION\n");
            break;
        case GL_OUT_OF_MEMORY: printf("GL_OUT_OF_MEMORY\n"); break;
        default: printf("<UNKNOWN GL ERROR>\n"); break;
    }
    return false;
}

//
// Clear window
//
void clear_buffers(const ym::vec4f& background) {
    assert(check_error());
    glClearColor(background[0], background[1], background[2], background[3]);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    assert(check_error());
}

//
// Enable/disable depth test
//
void enable_depth_test(bool enabled) {
    assert(check_error());
    if (enabled)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    assert(check_error());
}

//
// Enable/disable culling
//
void enable_culling(bool enabled) {
    assert(check_error());
    if (enabled)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
    assert(check_error());
}

//
// Enable/disable wireframe
//
void enable_wireframe(bool enabled) {
    assert(check_error());
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    assert(check_error());
}

//
// Enable/disable edges. Attempts to avoid z-fighting but the method is not
// robust.
//
void enable_edges(bool enabled, float tolerance) {
    assert(check_error());
    if (enabled) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDepthRange(0, tolerance);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDepthRange(0, 1);
    }
    assert(check_error());
}

//
// Enable/disable blending
//
void enable_blending(bool enabled) {
    assert(check_error());
    if (enabled) {
        glEnable(GL_BLEND);
    } else {
        glDisable(GL_BLEND);
    }
    assert(check_error());
}

//
// Set blending to over operator
//
void set_blend_over() {
    assert(check_error());
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    assert(check_error());
}

//
// Line width
//
void line_width(float w) {
    assert(check_error());
    glLineWidth(std::min(std::max(w, 0.0f), 1.0f));
    assert(check_error());
}

//
// Set viewport
//
void set_viewport(const ym::vec4i& v) {
    assert(check_error());
    glViewport(v.x, v.y, v.z, v.w);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void shade_image(
    uint tid, int win_w, int win_h, float ox, float oy, float zoom) {
    assert(check_error());
    shade_image(tid, win_w, win_h, ox, oy, zoom, ym::tonemap_type::none, 0, 1);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void shade_image(uint tid, int win_w, int win_h, float ox, float oy, float zoom,
    ym::tonemap_type tmtype, float exposure, float gamma) {
    static const std::string& header =
        R"(
            #version 330

            #define pi 3.14159265

            uniform vec2 offset;
            uniform vec2 win_size;
            uniform float zoom;

            uniform sampler2D img;

        )";

    static const std::string& vert =
        R"(
        layout(location = 0) in vec2 vert_texcoord;

        out vec2 texcoord;

        void main() {
            vec2 size = textureSize(img, 0).xy;
            texcoord = vert_texcoord.xy;
            vec2 pos = offset + size * vert_texcoord.xy * zoom;
            vec2 upos = 2 * pos / win_size - vec2(1,1);
            upos.y = - upos.y;
            gl_Position = vec4(upos.x, upos.y, 0, 1);
        }

        )";

    static const std::string frag_tonemap =
        R"(
        #define TONEMAP_LINEAR 0
        #define TONEMAP_SRGB 1
        #define TONEMAP_GAMMA 2
        #define TONEMAP_FILMIC 3

        struct Tonemap {
            int type;
            float exposure;
            float gamma;
        };
        uniform Tonemap tonemap;

        vec3 eval_filmic(vec3 x) {
            float a = 2.51f;
            float b = 0.03f;
            float c = 2.43f;
            float d = 0.59f;
            float e = 0.14f;
            return clamp((x*(a*x+b))/(x*(c*x+d)+e),0,1);
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tonemap.exposure);
            if(tonemap.type == TONEMAP_SRGB) {
                c = pow(c,vec3(1/2.2));
            } else if(tonemap.type == TONEMAP_GAMMA) {
                c = pow(c,vec3(1/tonemap.gamma));
            } else if(tonemap.type == TONEMAP_FILMIC) {
                c = eval_filmic(c);
            }
            return c;
        }

        )";

    static const std::string frag_main =
        R"(
        in vec2 texcoord;
        out vec4 color;

        void main() {
             vec4 c = texture(img,texcoord);
             c.xyz = eval_tonemap(c.xyz);
             color = c;
        }
        )";

    static uint fid, vid, vao;
    static uint prog = make_program(
        header + vert, header + frag_tonemap + frag_main, &vid, &fid, &vao);

    static uint tvbo_id = make_buffer(
        std::vector<ym::vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false, false);
    static uint evbo_id =
        make_buffer(std::vector<int>{0, 1, 2, 0, 2, 3}, true, false);

    assert(tid != 0);

    bind_vertex_arrays(vao);
    bind_program(prog);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tid);

    set_uniform(prog, "zoom", zoom);
    set_uniform(prog, "win_size", ym::vec2f{(float)win_w, (float)win_h});
    set_uniform(prog, "offset", ym::vec2f{ox, oy});
    set_uniform(prog, "tonemap.type", (int)tmtype);
    set_uniform(prog, "tonemap.exposure", exposure);
    set_uniform(prog, "tonemap.gamma", gamma);
    set_uniform_texture(prog, "img", tid, 0);

    set_vertattr(prog, "vert_texcoord", tvbo_id, ym::vec2f{0, 0});
    draw_elems(2, evbo_id, etype::triangle);

    bind_program(0);
    bind_vertex_arrays(0);

    glDisable(GL_BLEND);

    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void read_imagef(float* pixels, int w, int h, int nc) {
    assert(check_error());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
    assert(check_error());
}

//
// Implementation of make_texture.
//
uint _make_texture(int w, int h, int nc, const void* pixels, GLuint type,
    bool linear, bool mipmap, bool as_float, bool as_srgb) {
    assert(!as_srgb || !as_float);
    assert(check_error());
    int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
    int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
    int* formats =
        (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
    uint id;
    assert(check_error());
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0, formats_ub[nc - 1],
        type, pixels);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
    }
    glBindTexture(GL_TEXTURE_2D, 0);
    assert(check_error());
    return id;
}

//
// This is a public API. See above for documentation.
//
uint make_texture(int w, int h, int nc, const float* pixels, bool linear,
    bool mipmap, bool as_float) {
    return _make_texture(
        w, h, nc, pixels, GL_FLOAT, linear, mipmap, as_float, false);
}

//
// This is a public API. See above for documentation.
//
uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap, bool as_srgb) {
    return _make_texture(
        w, h, nc, pixels, GL_UNSIGNED_BYTE, linear, mipmap, false, as_srgb);
}

//
// Implementation of update_texture.
//
void _update_texture(uint id, int w, int h, int nc, const void* pixels,
    GLuint type, bool mipmap) {
    assert(check_error());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(
        GL_TEXTURE_2D, 0, 0, 0, w, h, formats[nc - 1], type, pixels);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void update_texture(
    uint id, int w, int h, int nc, const float* pixels, bool mipmap) {
    assert(check_error());
    _update_texture(id, w, h, nc, pixels, GL_FLOAT, mipmap);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void update_texture(
    uint id, int w, int h, int nc, const unsigned char* pixels, bool mipmap) {
    assert(check_error());
    _update_texture(id, w, h, nc, pixels, GL_UNSIGNED_BYTE, mipmap);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void clear_texture(uint* tid) {
    assert(check_error());
    glDeleteTextures(1, tid);
    *tid = 0;
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
uint make_buffer_raw(
    int size, const void* values, buffer_type btype, bool dynamic) {
    assert(check_error());
    static auto type_map = std::map<buffer_type, int>{
        {buffer_type::vertex, GL_ARRAY_BUFFER},
        {buffer_type::element, GL_ELEMENT_ARRAY_BUFFER},
        {buffer_type::uniform_block, GL_UNIFORM_BUFFER},
    };
    auto target = type_map.at(btype);
    auto bid = (GLuint)0;
    glGenBuffers(1, &bid);
    glBindBuffer(target, bid);
    glBufferData(
        target, size, values, (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    glBindBuffer(target, 0);
    return bid;
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void update_buffer_raw(
    uint bid, int size, const void* values, buffer_type btype, bool dynamic) {
    assert(check_error());
    static auto type_map = std::map<buffer_type, int>{
        {buffer_type::vertex, GL_ARRAY_BUFFER},
        {buffer_type::element, GL_ELEMENT_ARRAY_BUFFER},
        {buffer_type::uniform_block, GL_UNIFORM_BUFFER},
    };
    auto target = type_map.at(btype);
    glBindBuffer(target, bid);
    glBufferSubData(target, 0, size, values);
    glBindBuffer(target, 0);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void clear_buffer(uint* bid) {
    assert(check_error());
    glDeleteBuffers(1, bid);
    *bid = 0;
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
uint make_vertex_arrays() {
    assert(check_error());
    auto aid = (uint)0;
    glGenVertexArrays(1, &aid);
    assert(check_error());
    return aid;
}

//
// This is a public API. See above for documentation.
//
void clear_vertex_arrays(uint* aid) {
    assert(check_error());
    glDeleteVertexArrays(1, aid);
    *aid = 0;
    assert(check_error());
}

//
// Binds a vertex array object
//
void bind_vertex_arrays(uint aid) { glBindVertexArray(aid); }

//
// This is a public API. See above for documentation.
//
uint make_program(const std::string& vertex, const std::string& fragment,
    uint* vid, uint* fid, uint* vao) {
    assert(check_error());
    uint gl_vao = 0;
    glGenVertexArrays(1, &gl_vao);
    glBindVertexArray(gl_vao);
    assert(check_error());

    int errflags[2];
    char errbuf[10000];
    int gl_shader[2] = {0, 0};
    std::string code[2] = {vertex, fragment};
    int type[2] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER};
    for (int i = 0; i < 2; i++) {
        gl_shader[i] = glCreateShader(type[i]);
        const char* code_str = code[i].c_str();
        glShaderSource(gl_shader[i], 1, &code_str, NULL);
        glCompileShader(gl_shader[i]);
        glGetShaderiv(gl_shader[i], GL_COMPILE_STATUS, &errflags[i]);
        if (!errflags[i]) {
            glGetShaderInfoLog(gl_shader[i], 10000, 0, errbuf);
            printf("shader not compiled\n\n%s\n\n", errbuf);
            return 0;
        }
        assert(glGetError() == GL_NO_ERROR);
    }

    // create program
    int gl_program = glCreateProgram();
    for (int i = 0; i < 2; i++) glAttachShader(gl_program, gl_shader[i]);
    glLinkProgram(gl_program);
    glValidateProgram(gl_program);
    glGetProgramiv(gl_program, GL_LINK_STATUS, &errflags[0]);
    glGetProgramiv(gl_program, GL_VALIDATE_STATUS, &errflags[1]);
    if (!errflags[0] || !errflags[1]) {
        glGetProgramInfoLog(gl_program, 10000, 0, errbuf);
        printf("program not linked\n\n%s\n\n", errbuf);
        return 0;
    }
    assert(check_error());

    glBindVertexArray(0);
    assert(check_error());

    // returns
    if (vid) *vid = gl_shader[0];
    if (fid) *fid = gl_shader[1];
    if (vao) *vao = gl_vao;
    assert(check_error());
    return gl_program;
}

//
// This is a public API. See above for documentation.
//
void clear_program(uint* pid, uint* vid, uint* fid, uint* vao) {
    assert(check_error());
    if (vid) {
        glDetachShader(*pid, *vid);
        glDeleteShader(*vid);
        *vid = 0;
    }
    if (fid) {
        glDetachShader(*pid, *fid);
        glDeleteShader(*fid);
        *fid = 0;
    }
    if (pid) {
        glDeleteProgram(*pid);
        *pid = 0;
    }
    if (vao) {
        glDeleteVertexArrays(1, vao);
        *vao = false;
    }
    assert(check_error());
}

//
// Binds a program object
//
void bind_program(uint pid) { glUseProgram(pid); }

//
// Get uniform location (simple GL wrapper that avoids GL includes)
//
int get_uniform_location(uint pid, const std::string& name) {
    assert(check_error());
    return glGetUniformLocation(pid, name.c_str());
    assert(check_error());
}

//
// Get attrib location (simple GL wrapper that avoids GL includes)
//
int get_attrib_location(uint pid, const std::string& name) {
    assert(check_error());
    return glGetAttribLocation(pid, name.c_str());
    assert(check_error());
}

//
// Get the names of all uniforms
//
std::vector<std::pair<std::string, int>> get_uniforms_names(uint pid) {
    auto num = 0;
    assert(check_error());
    glGetProgramiv(pid, GL_ACTIVE_UNIFORMS, &num);
    assert(check_error());
    auto names = std::vector<std::pair<std::string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0, length = 0;
        GLenum type;
        glGetActiveUniform(pid, i, 4096, &length, &size, &type, name);
        if (length > 3 && name[length - 1] == ']' && name[length - 2] == '0' &&
            name[length - 3] == '[')
            name[length - 3] = 0;
        auto loc = glGetUniformLocation(pid, name);
        if (loc < 0) continue;
        names.push_back({name, loc});
        assert(check_error());
    }
    return names;
}

//
// Get the names of all uniforms
//
std::vector<std::pair<std::string, int>> get_uniform_buffers_names(uint pid) {
    auto num = 0;
    assert(check_error());
    glGetProgramiv(pid, GL_ACTIVE_UNIFORM_BLOCKS, &num);
    assert(check_error());
    auto names = std::vector<std::pair<std::string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        glGetActiveUniformBlockName(pid, i, 4096, nullptr, name);
        auto loc = glGetUniformBlockIndex(pid, name);
        names.push_back({name, loc});
        assert(check_error());
    }
    return names;
}

//
// Get the names of all attributes
//
std::vector<std::pair<std::string, int>> get_attributes_names(uint pid) {
    auto num = 0;
    assert(check_error());
    glGetProgramiv(pid, GL_ACTIVE_ATTRIBUTES, &num);
    assert(check_error());
    auto names = std::vector<std::pair<std::string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0;
        GLenum type;
        glGetActiveAttrib(pid, i, 4096, nullptr, &size, &type, name);
        auto loc = glGetAttribLocation(pid, name);
        if (loc < 0) continue;
        names.push_back({name, loc});
        assert(check_error());
    }
    return names;
}

//
// This is a public API. See above for documentation.
//
bool set_uniform_raw(uint prog, int pos, const int* val, int nc, int count) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    if (pos < 0) return false;
    switch (nc) {
        case 1: glUniform1iv(pos, count, val); break;
        case 2: glUniform2iv(pos, count, val); break;
        case 3: glUniform3iv(pos, count, val); break;
        case 4: glUniform4iv(pos, count, val); break;
        default: assert(false);
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_uniform_raw(uint prog, int pos, const float* val, int nc, int count) {
    assert((nc >= 1 && nc <= 4) || (nc == 16) || (nc == 12));
    assert(check_error());
    if (pos < 0) return false;
    switch (nc) {
        case 1: glUniform1fv(pos, count, val); break;
        case 2: glUniform2fv(pos, count, val); break;
        case 3: glUniform3fv(pos, count, val); break;
        case 4: glUniform4fv(pos, count, val); break;
        case 12: glUniformMatrix4x3fv(pos, count, false, val); break;
        case 16: glUniformMatrix4fv(pos, count, false, val); break;
        default: assert(false); return 0;
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_uniform_block(uint prog, int pos, int bid, uint bunit) {
    assert(check_error());
    if (pos < 0) return false;
    if (bid < 0) return false;
    glBindBufferBase(GL_UNIFORM_BUFFER, bunit, bid);
    glUniformBlockBinding(prog, pos, bunit);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_uniform_texture(
    uint prog, int pos, const texture_info& tinfo, uint tunit) {
    static const auto wrap_mode_map =
        std::map<texture_wrap, uint>{{texture_wrap::repeat, GL_REPEAT},
            {texture_wrap::clamp, GL_CLAMP_TO_EDGE},
            {texture_wrap::mirror, GL_MIRRORED_REPEAT}};
    static const auto filter_mode_map =
        std::map<texture_filter, uint>{{texture_filter::nearest, GL_NEAREST},
            {texture_filter::linear, GL_LINEAR},
            {texture_filter::nearest_mipmap_nearest, GL_NEAREST_MIPMAP_NEAREST},
            {texture_filter::linear_mipmap_nearest, GL_LINEAR_MIPMAP_NEAREST},
            {texture_filter::nearest_mipmap_linear, GL_NEAREST_MIPMAP_LINEAR},
            {texture_filter::linear_mipmap_linear, GL_LINEAR_MIPMAP_LINEAR}};

    assert(check_error());
    if (pos < 0) return false;
    if (tinfo.txt_id > 0) {
        glActiveTexture(GL_TEXTURE0 + tunit);
        glBindTexture(GL_TEXTURE_2D, tinfo.txt_id);
        if (tinfo.wrap_s != texture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                wrap_mode_map.at(tinfo.wrap_s));
        if (tinfo.wrap_t != texture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                wrap_mode_map.at(tinfo.wrap_t));
        if (tinfo.filter_min != texture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                filter_mode_map.at(tinfo.filter_min));
        if (tinfo.filter_mag != texture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                filter_mode_map.at(tinfo.filter_mag));
        glUniform1i(pos, tunit);
    } else {
        glActiveTexture(GL_TEXTURE0 + tunit);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(pos, tunit);
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_ptr(
    uint prog, const std::string& var, const float* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    if (value) {
        glEnableVertexAttribArray(pos);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, value);
    } else {
        glDisableVertexAttribArray(pos);
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_buffer(uint prog, const std::string& var, uint bid, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    if (bid) {
        glEnableVertexAttribArray(pos);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        glDisableVertexAttribArray(pos);
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattri_buffer(uint prog, const std::string& var, uint bid, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    if (bid) {
        glEnableVertexAttribArray(pos);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        // glVertexAttribPointer(pos, nc, GL_INT, false, 0, 0);
        glVertexAttribIPointer(pos, nc, GL_INT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        glDisableVertexAttribArray(pos);
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_val_raw(uint prog, int pos, const float* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    if (pos < 0) return false;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttrib1fv(pos, value); break;
        case 2: glVertexAttrib2fv(pos, value); break;
        case 3: glVertexAttrib3fv(pos, value); break;
        case 4: glVertexAttrib4fv(pos, value); break;
        default: assert(false); break;
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_val_raw(uint prog, int pos, const int* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    if (pos < 0) return false;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttribI1iv(pos, value); break;
        case 2: glVertexAttribI2iv(pos, value); break;
        case 3: glVertexAttribI3iv(pos, value); break;
        case 4: glVertexAttribI4iv(pos, value); break;
        default: assert(false); break;
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_raw(uint prog, int pos, uint bid, int nc, const float* def) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    if (pos < 0) return false;
    if (bid) {
        glEnableVertexAttribArray(pos);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        glDisableVertexAttribArray(pos);
        if (def) {
            switch (nc) {
                case 1: glVertexAttrib1fv(pos, def); break;
                case 2: glVertexAttrib2fv(pos, def); break;
                case 3: glVertexAttrib3fv(pos, def); break;
                case 4: glVertexAttrib4fv(pos, def); break;
                default: assert(false); break;
            }
        }
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool set_vertattr_raw(uint prog, int pos, uint bid, int nc, const int* def) {
    assert(nc >= 1 && nc <= 4);
    assert(check_error());
    if (pos < 0) return false;
    if (bid) {
        glEnableVertexAttribArray(pos);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        // glVertexAttribPointer(pos, nc, GL_INT, false, 0, 0);
        glVertexAttribIPointer(pos, nc, GL_INT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        glDisableVertexAttribArray(pos);
        if (def) {
            switch (nc) {
                case 1: glVertexAttribI1iv(pos, def); break;
                case 2: glVertexAttribI2iv(pos, def); break;
                case 3: glVertexAttribI3iv(pos, def); break;
                case 4: glVertexAttribI4iv(pos, def); break;
                default: assert(false); break;
            }
        }
    }
    assert(check_error());
    return true;
}

//
// This is a public API. See above for documentation.
//
bool draw_elems(int nelems, uint bid, etype type) {
    if (!nelems) return true;
    assert(bid);
    assert(check_error());
    int mode = 0, nc = 0;
    switch (type) {
        case etype::point:
            mode = GL_POINTS;
            nc = 1;
            break;
        case etype::line:
            mode = GL_LINES;
            nc = 2;
            break;
        case etype::triangle:
            mode = GL_TRIANGLES;
            nc = 3;
            break;
        case etype::quad:
            mode = GL_QUADS;
            nc = 4;
            break;
        default: assert(false);
    };
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(mode, nelems * nc, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(check_error());
    return true;
}

namespace stdshader {

//
// This is a public API. See above for documentation.
//
// Filmic tone mapping from
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
//
uint make_program(uint* vao) {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    static const std::string& vert_header =
        R"(
        #version 330

    )";

    static const std::string& vert_skinning =
        R"(
            #define SKIN_NONE 0
            #define SKIN_STD 1
            #define SKIN_GLTF 2

    uniform int skin_type = 0;
    uniform mat4 skin_xforms[32];
    layout(location = 6) in vec4 vert_skin_weights;            // vertex skinning weights
    layout(location = 7) in ivec4 vert_skin_joints;            // vertex skinning joints (in mesh coordinate frame)

    vec3 transform_point(mat4 m, vec3 p) {
        vec4 p4 = m * vec4(p,1);
        return p4.xyz / p4.w;
    }

    vec3 transform_normal(mat4 m, vec3 p) {
        vec4 p4 = m * vec4(p,0);
        return p4.xyz;
    }

    void apply_skin(inout vec3 pos, inout vec3 norm) {
        if(skin_type == 0) {
            return;
        } else if(skin_type == SKIN_STD) {
            vec4 w = vert_skin_weights;
            ivec4 j = ivec4(vert_skin_joints);
            pos = transform_point( skin_xforms[j.x], pos ) * w.x +
                  transform_point( skin_xforms[j.y], pos ) * w.y +
                  transform_point( skin_xforms[j.z], pos ) * w.z +
                  transform_point( skin_xforms[j.w], pos ) * w.w;
            norm = normalize(
                   transform_normal( skin_xforms[j.x], norm ) * w.x +
                   transform_normal( skin_xforms[j.y], norm ) * w.y +
                   transform_normal( skin_xforms[j.z], norm ) * w.z +
                   transform_normal( skin_xforms[j.w], norm ) * w.w);
         } else if(skin_type == SKIN_GLTF) {
             vec4 w = vert_skin_weights;
             ivec4 j = ivec4(vert_skin_joints);
             mat4 xf = skin_xforms[j.x] * w.x + skin_xforms[j.y] * w.y +
                       skin_xforms[j.z] * w.z + skin_xforms[j.w] * w.w;
            pos = transform_point(xf, pos);
            norm = normalize(transform_normal(xf, norm));
         }
    }
    )";

    static const std::string& vert_main =
        R"(
    layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
    layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
    layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
    layout(location = 3) in vec4 vert_color;          // vertex color
    layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

    uniform mat4 shape_xform;           // shape transform

    struct Camera {
        mat4 xform;          // camera xform
        mat4 xform_inv;      // inverse of the camera frame (as a matrix)
        mat4 proj;           // camera projection
    };
    uniform Camera camera;      // camera data

    out vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
    out vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
    out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
    out vec4 color;                 // [to fragment shader] vertex color
    out vec4 tangsp;                // [to fragment shader] vertex tangent space

    // main function
    void main() {
        // copy values
        pos = vert_pos;
        norm = vert_norm;
        tangsp = vert_tangsp;

        // world projection
        pos = (shape_xform * vec4(pos,1)).xyz;
        norm = (shape_xform * vec4(norm,0)).xyz;
        tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

        // skinning
        apply_skin(pos, norm);

        // copy other vertex properties
        texcoord = vert_texcoord;
        color = vert_color;

        // clip
        gl_Position = camera.proj * camera.xform_inv * vec4(pos,1);
    }
    )";

    static const std::string frag_header =
        R"(
        #version 330

        #define pi 3.14159265

        )";

    static const std::string frag_tonemap =
        R"(
        #define TONEMAP_LINEAR 0
        #define TONEMAP_SRGB 1
        #define TONEMAP_GAMMA 2
        #define TONEMAP_FILMIC 3

        struct Tonemap {
            int type;       // tonemap type (TM_...)
            float exposure; // image exposure
            float gamma;    // image gamma
        };
        uniform Tonemap tonemap;

        vec3 eval_filmic(vec3 x) {
            float a = 2.51f;
            float b = 0.03f;
            float c = 2.43f;
            float d = 0.59f;
            float e = 0.14f;
            return clamp((x*(a*x+b))/(x*(c*x+d)+e),0,1);
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tonemap.exposure);
            if(tonemap.type == TONEMAP_SRGB) {
                c = pow(c,vec3(1/2.2));
            } else if(tonemap.type == TONEMAP_GAMMA) {
                c = pow(c,vec3(1/tonemap.gamma));
            } else if(tonemap.type == TONEMAP_FILMIC) {
                c = eval_filmic(c);
            }
            return c;
        }

        )";

    static const std::string frag_lighting =
        R"(
        struct Lighting {
            bool eyelight;        // eyelight shading
            vec3 amb;             // ambient light
            int lnum;              // number of lights
            int ltype[16];         // light type (0 -> point, 1 -> directional)
            vec3 lpos[16];         // light positions
            vec3 lke[16];          // light intensities
        };
        uniform Lighting lighting;

        void eval_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(lighting.ltype[lid] == 0) {
                // compute point light color at pos
                cl = lighting.lke[lid] / pow(length(lighting.lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lighting.lpos[lid]-pos);
            }
            else if(lighting.ltype[lid] == 1) {
                // compute light color
                cl = lighting.lke[lid];
                // compute light direction
                wi = normalize(lighting.lpos[lid]);
            }
        }

        )";

    static const std::string frag_brdf =
        R"(
        #define BRDF_NONE 0
        #define BRDF_POINT 1
        #define BRDF_KAJIYA_KAY 2
        #define BRDF_GGX 3
        #define BRDF_PHONG 4

        struct Brdf {
            int type;
            vec3 ke;
            vec3 kd;
            vec3 ks;
            float rs;
            float op;
            bool cutout;
        };

        vec3 brdfcos(Brdf brdf, vec3 n, vec3 wi, vec3 wo) {
            if(brdf.type == BRDF_NONE) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(brdf.rs*brdf.rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(brdf.type == BRDF_POINT) {
                return ((1+dot(wo,wi))/2) * brdf.kd/pi;
            } else if(brdf.type == BRDF_KAJIYA_KAY) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * brdf.kd / pi;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pi)) * pow(si,ns);
                vec3 spec = si * brdf.ks * d / (4*si*so);
                return diff+spec;
            } else if(brdf.type == BRDF_GGX || brdf.type == BRDF_PHONG) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * brdf.kd / pi;
                if(ndh<=0) return diff;
                if(brdf.type == BRDF_PHONG) {
                    float d = ((2+ns)/(2*pi)) * pow(ndh,ns);
                    vec3 spec = ndi * brdf.ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = brdf.rs * brdf.rs;
                    float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * brdf.ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        )";

    static const std::string frag_material =
        R"(
        #define MATERIAL_EMISSION_ONLY 0
        #define MATERIAL_GENERIC 1
        #define MATERIAL_GLTF_METALLIC_ROUGHNESS 2
        #define MATERIAL_GLTF_SPECULAR_GLOSSINESS 3

        #define ELEMENT_POINT 1
        #define ELEMENT_LINE 2
        #define ELEMENT_TRIANGLE 3

        struct Material {
            int mtype;         // material type
            int etype;         // element type
            vec3 ke;           // material ke
            vec3 kd;           // material kd
            vec3 ks;           // material ks
            float rs;          // material rs
            float op;          // material op

            bool txt_ke_on;    // material ke texture on
            sampler2D txt_ke;  // material ke texture
            bool txt_kd_on;    // material kd texture on
            sampler2D txt_kd;  // material kd texture
            bool txt_ks_on;    // material ks texture on
            sampler2D txt_ks;  // material ks texture
            bool txt_rs_on;    // material rs texture on
            sampler2D txt_rs;  // material rs texture

            bool txt_norm_on;    // material norm texture on
            sampler2D txt_norm;  // material norm texture
            sampler2D txt_norm_scale;  // material norm scale

            bool txt_occ_on;    // material occ texture on
            sampler2D txt_occ;  // material occ texture
            sampler2D txt_occ_scale;  // material occ scale

            bool use_phong;       // material use phong
            bool double_sided;    // material double sided
            bool alpha_cutout;    // material alpha cutout
        };
        uniform Material material;

        void eval_material(vec2 texcoord, vec4 color, out int type, out vec3 ke,
                           out vec3 kd, out vec3 ks, out float rs, out float op, out bool cutout) {
            ke = color.xyz * material.ke;
            kd = color.xyz * material.kd;
            ks = color.xyz * material.ks;
            rs = material.rs;
            op = color.w * material.op;

            vec3 ke_txt = (material.txt_ke_on) ? texture(material.txt_ke,texcoord).xyz : vec3(1);
            vec4 kd_txt = (material.txt_kd_on) ? texture(material.txt_kd,texcoord) : vec4(1);
            vec4 ks_txt = (material.txt_ks_on) ? texture(material.txt_ks,texcoord) : vec4(1);
            float rs_txt = (material.txt_rs_on) ? texture(material.txt_rs,texcoord).x : 1;

            // scale common values
            ke *= ke_txt;

            // get material color from textures and adjust values
            if(material.mtype == MATERIAL_EMISSION_ONLY) {
                type = BRDF_NONE;
            } else if(material.mtype == MATERIAL_GENERIC) {
                type = material.etype;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt;
            } else if(material.mtype == MATERIAL_GLTF_METALLIC_ROUGHNESS) {
                type = material.etype;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(material.mtype == MATERIAL_GLTF_SPECULAR_GLOSSINESS) {
                type = material.etype;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= ks_txt.w;
                rs = (1 - rs) * (1 - rs);
                op *= kd_txt.w;
            }

            cutout = material.alpha_cutout && op == 0;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!material.txt_norm_on) return norm;
            vec3 tangu = normalize(tangsp.xyz);
            vec3 tangv = normalize(cross(tangu, norm));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 txt = 2 * pow(texture(material.txt_norm,texcoord).xyz, vec3(1/2.2)) - 1;
            return normalize( tangu * txt.x + tangv * txt.y + norm * txt.z );
        }

        )";

    static const std::string frag_main =
        R"(
        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        struct Camera {
            mat4 xform;          // camera xform
            mat4 xform_inv;      // inverse of the camera frame (as a matrix)
            mat4 proj;           // camera projection
        };
        uniform Camera camera;      // camera data

        uniform vec4 highlight;   // highlighted color

        out vec4 frag_color;        // eyelight shading

        // main
        void main() {
            // view vector
            vec3 wo = normalize( (camera.xform*vec4(0,0,0,1)).xyz - pos );

            // re-normalize normals
            vec3 n = normalize(norm);

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(material.double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            Brdf brdf;
            eval_material(texcoord, color, brdf.type, brdf.ke, brdf.kd, brdf.ks, brdf.rs, brdf.op, brdf.cutout);

            // exit if needed
            if(brdf.cutout) discard;

            // emission
            vec3 c = brdf.ke;

            // check early exit
            if(brdf.kd != vec3(0,0,0) || brdf.ks != vec3(0,0,0)) {
                // eyelight shading
                if(lighting.eyelight) {
                    vec3 wi = wo;
                    c += pi * brdfcos(brdf,n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lighting.amb * brdf.kd;
                    // foreach light
                    for(int lid = 0; lid < lighting.lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        eval_light(lid, pos, cl, wi);
                        c += cl * brdfcos(brdf,n,wi,wo);
                    }
                }
            }

            // final color correction
            c = eval_tonemap(c);

            // highlighting
            if(highlight.w > 0) {
                if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
                    c = highlight.xyz * highlight.w + c * (1-highlight.w);
            }

            // output final color by setting gl_FragColor
            frag_color = vec4(c,brdf.op);
        }
    )";

    assert(check_error());
    uint vid, fid;
    return yglu::make_program(vert_header + vert_skinning + vert_main,
        frag_header + frag_tonemap + frag_lighting + frag_brdf + frag_material +
            frag_main,
        &vid, &fid, vao);
    assert(check_error());
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif
}

//
// This is a public API. See above for documentation.
//
void begin_frame(uint prog, uint vao, bool shade_eyelight, float img_exposure,
    ym::tonemap_type img_tonemap, float img_gamma,
    const ym::mat4f& camera_xform, const ym::mat4f& camera_xform_inv,
    const ym::mat4f& camera_proj) {
    static auto eyelight_id = get_uniform_location(prog, "lighting.eyelight");
    static auto exposure_id = get_uniform_location(prog, "tonemap.exposure");
    static auto gamma_id = get_uniform_location(prog, "tonemap.gamma");
    static auto type_id = get_uniform_location(prog, "tonemap.type");
    static auto xform_id = get_uniform_location(prog, "camera.xform");
    static auto xform_inv_id = get_uniform_location(prog, "camera.xform_inv");
    static auto proj_id = get_uniform_location(prog, "camera.proj");
    assert(check_error());
    bind_vertex_arrays(vao);
    bind_program(prog);
    set_uniform(prog, eyelight_id, shade_eyelight);
    set_uniform(prog, exposure_id, img_exposure);
    set_uniform(prog, gamma_id, img_gamma);
    set_uniform(prog, type_id, (int)img_tonemap);
    set_uniform(prog, xform_id, camera_xform);
    set_uniform(prog, xform_inv_id, camera_xform_inv);
    set_uniform(prog, proj_id, camera_proj);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void end_frame() {
    assert(check_error());
    glBindVertexArray(0);
    glUseProgram(0);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void set_lights(uint prog, const ym::vec3f& amb, int num, ym::vec3f* pos,
    ym::vec3f* ke, ltype* type) {
    static auto amb_id = get_uniform_location(prog, "lighting.amb");
    static auto lnum_id = get_uniform_location(prog, "lighting.lnum");
    static auto lpos_id = get_uniform_location(prog, "lighting.lpos");
    static auto lke_id = get_uniform_location(prog, "lighting.lke");
    static auto ltype_id = get_uniform_location(prog, "lighting.ltype");
    assert(check_error());
    set_uniform(prog, amb_id, amb);
    set_uniform(prog, lnum_id, num);
    set_uniform(prog, lpos_id, pos, num);
    set_uniform(prog, lke_id, ke, num);
    set_uniform(prog, ltype_id, (int*)type, num);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void begin_shape(uint prog, const ym::mat4f& xform) {
    static auto xform_id = get_uniform_location(prog, "shape_xform");
    assert(check_error());
    set_uniform(prog, xform_id, xform);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
void end_shape() {
    assert(check_error());
    for (int i = 0; i < 16; i++) glDisableVertexAttribArray(i);
    assert(check_error());
}

//
// Set the object as highlighted.
//
void set_highlight(uint prog, const ym::vec4f& highlight) {
    static auto highlight_id = get_uniform_location(prog, "highlight");
    set_uniform(prog, highlight_id, highlight);
}

//
// Internal representation
//
inline void set_material_generic(uint prog, int mtype, int etype,
    const ym::vec3f& ke, const ym::vec3f& kd, const ym::vec3f& ks, float rs,
    float op, const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& rs_txt,
    const texture_info& norm_txt, const texture_info& occ_txt, bool use_phong,
    bool double_sided, bool alpha_cutout) {
    static auto mtype_id = get_uniform_location(prog, "material.mtype");
    static auto etype_id = get_uniform_location(prog, "material.etype");
    static auto ke_id = get_uniform_location(prog, "material.ke");
    static auto kd_id = get_uniform_location(prog, "material.kd");
    static auto ks_id = get_uniform_location(prog, "material.ks");
    static auto rs_id = get_uniform_location(prog, "material.rs");
    static auto op_id = get_uniform_location(prog, "material.op");
    static auto ke_txt_id = get_uniform_location(prog, "material.txt_ke");
    static auto ke_txt_on_id = get_uniform_location(prog, "material.txt_ke_on");
    static auto kd_txt_id = get_uniform_location(prog, "material.txt_kd");
    static auto kd_txt_on_id = get_uniform_location(prog, "material.txt_kd_on");
    static auto ks_txt_id = get_uniform_location(prog, "material.txt_ks");
    static auto ks_txt_on_id = get_uniform_location(prog, "material.txt_ks_on");
    static auto rs_txt_id = get_uniform_location(prog, "material.txt_rs");
    static auto rs_txt_on_id = get_uniform_location(prog, "material.txt_rs_on");
    static auto norm_txt_id = get_uniform_location(prog, "material.txt_norm");
    static auto norm_txt_on_id =
        get_uniform_location(prog, "material.txt_norm_on");
    static auto occ_txt_id = get_uniform_location(prog, "material.txt_occ");
    static auto occ_txt_on_id =
        get_uniform_location(prog, "material.txt_occ_on");
    static auto norm_scale_id =
        get_uniform_location(prog, "material.norm_scale");
    static auto occ_scale_id = get_uniform_location(prog, "material.occ_scale");
    static auto use_phong_id = get_uniform_location(prog, "material.use_phong");
    static auto double_sided_id =
        get_uniform_location(prog, "material.double_sided");
    static auto alpha_cutout_id =
        get_uniform_location(prog, "material.alpha_cutout");

    assert(check_error());
    set_uniform(prog, mtype_id, mtype);
    set_uniform(prog, etype_id, (int)etype);
    set_uniform(prog, ke_id, ke);
    set_uniform(prog, kd_id, kd);
    set_uniform(prog, ks_id, ks);
    set_uniform(prog, rs_id, rs);
    set_uniform(prog, op_id, op);
    set_uniform_texture(prog, ke_txt_id, ke_txt_on_id, ke_txt, 0);
    set_uniform_texture(prog, kd_txt_id, kd_txt_on_id, kd_txt, 1);
    set_uniform_texture(prog, ks_txt_id, ks_txt_on_id, ks_txt, 2);
    set_uniform_texture(prog, rs_txt_id, rs_txt_on_id, rs_txt, 3);
    set_uniform_texture(prog, norm_txt_id, norm_txt_on_id, norm_txt, 4);
    set_uniform_texture(prog, occ_txt_id, occ_txt_on_id, occ_txt, 5);
    set_uniform(prog, norm_scale_id, norm_txt.scale);
    set_uniform(prog, occ_scale_id, occ_txt.scale);
    set_uniform(prog, use_phong_id, use_phong);
    set_uniform(prog, double_sided_id, double_sided);
    set_uniform(prog, alpha_cutout_id, alpha_cutout);
    assert(check_error());
}

//
// Set material values for emission only (constant color).
// Indicates textures ids with the correspoinding XXX_txt variables.
//
void set_material_emission_only(uint prog, etype etype, const ym::vec3f& ke,
    float op, const texture_info& ke_txt, bool double_sided,
    bool alpha_cutout) {
    set_material_generic(prog, 0, (int)etype, ke, {0, 0, 0}, {0, 0, 0}, 1, op,
        ke_txt, 0, 0, 0, 0, 0, false, double_sided, alpha_cutout);
}

//
// Set material values with emission ke, diffuse kd, specular ks and
// specular roughness rs, opacity op. Indicates textures ids with the
// correspoinding
// XXX_txt variables. Uses GGX by default, but can switch to Phong is needed.
//
void set_material_generic(uint prog, etype etype, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op,
    const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& rs_txt,
    const texture_info& norm_txt, const texture_info& occ_txt, bool use_phong,
    bool double_sided, bool alpha_cutout) {
    set_material_generic(prog, 1, (int)etype, ke, kd, ks, rs, op, ke_txt,
        kd_txt, ks_txt, rs_txt, norm_txt, occ_txt, use_phong, double_sided,
        alpha_cutout);
}

//
// Set material values for glTF specular-roughness PBR shader,
// with emission ke, base color kb, opacity op, metallicity km and
// specular roughness rs. Uses basecolor-opacity texture kb_txt and
// metallic-roughness texture km_txt. Uses GGX but can be switched to Phong.
//
void set_material_gltf_metallic_roughness(uint prog, etype etype,
    const ym::vec3f& ke, const ym::vec3f& kb, float km, float rs, float op,
    const texture_info& ke_txt, const texture_info& kb_txt,
    const texture_info& km_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout) {
    set_material_generic(prog, 2, (int)etype, ke, kb, {km, km, km}, rs, op,
        ke_txt, kb_txt, km_txt, 0, norm_txt, occ_txt, use_phong, double_sided,
        alpha_cutout);
}

//
// Set material values for glTF specular-roughness PBR shader,
// with emission ke, diffuse color kd, opacity op, specular ks and
// specular glossiness rs. Uses diffuse-opacity texture kd_txt and
// specular-glpossiness texture ks_txt. Uses GGX but can be switched to Phong.
//
void set_material_gltf_specular_glossiness(uint prog, etype etype,
    const ym::vec3f& ke, const ym::vec3f& kd, const ym::vec3f& ks, float rs,
    float op, const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout) {
    set_material_generic(prog, 3, (int)etype, ke, kd, ks, rs, op, ke_txt,
        kd_txt, ks_txt, 0, norm_txt, occ_txt, use_phong, double_sided,
        alpha_cutout);
}

//
// This is a public API. See above for documentation.
//
float specular_exponent_to_roughness(float n) { return std::sqrt(2 / (n + 2)); }

//
// This is a public API. See above for documentation.
//
void set_vert(
    uint prog, uint pos, uint norm, uint texcoord, uint color, uint tangsp) {
    static auto pos_id = get_attrib_location(prog, "vert_pos");
    static auto norm_id = get_attrib_location(prog, "vert_norm");
    static auto texcoord_id = get_attrib_location(prog, "vert_texcoord");
    static auto color_id = get_attrib_location(prog, "vert_color");
    static auto tangsp_id = get_attrib_location(prog, "vert_tangsp");
    assert(check_error());
    set_vertattr(prog, pos_id, pos, ym::zero3f);
    set_vertattr(prog, norm_id, norm, ym::zero3f);
    set_vertattr(prog, texcoord_id, texcoord, ym::zero2f);
    set_vertattr(prog, color_id, color, ym::one4f);
    set_vertattr(prog, tangsp_id, tangsp, ym::zero4f);
    assert(check_error());
}

//
// Set vertex data with buffers for skinning.
//
void set_vert_skinning(uint prog, uint weights, uint joints, int nxforms,
    const ym::mat4f* xforms) {
    static auto type_id = get_uniform_location(prog, "skin_type");
    static auto xforms_id = get_uniform_location(prog, "skin_xforms");
    static auto weights_id = get_attrib_location(prog, "vert_skin_weights");
    static auto joints_id = get_attrib_location(prog, "vert_skin_joints");
    int type = 1;
    set_uniform(prog, type_id, type);
    set_uniform(prog, xforms_id, xforms, std::min(nxforms, 32));
    set_vertattr(prog, weights_id, weights, ym::zero4f);
    set_vertattr(prog, joints_id, joints, ym::zero4i);
}

//
// Set vertex data with buffers for skinning.
//
void set_vert_gltf_skinning(uint prog, uint weights, uint joints, int nxforms,
    const ym::mat4f* xforms) {
    static auto type_id = get_uniform_location(prog, "skin_type");
    static auto xforms_id = get_uniform_location(prog, "skin_xforms");
    static auto weights_id = get_attrib_location(prog, "vert_skin_weights");
    static auto joints_id = get_attrib_location(prog, "vert_skin_joints");
    int type = 2;
    set_uniform(prog, type_id, type);
    set_uniform(prog, xforms_id, xforms, std::min(nxforms, 32));
    set_vertattr(prog, weights_id, weights, ym::zero4f);
    set_vertattr(prog, joints_id, joints, ym::zero4i);
}

//
// Disables vertex skinning.
//
void set_vert_skinning_off(uint prog) {
    static auto type_id = get_uniform_location(prog, "skin_type");
    // static auto xforms_id = get_uniform_location(prog, "skin_xforms");
    static auto weights_id = get_attrib_location(prog, "vert_skin_weights");
    static auto joints_id = get_attrib_location(prog, "vert_skin_joints");
    int type = 0;
    set_uniform(prog, type_id, type);
    set_vertattr(prog, weights_id, 0, ym::zero4f);
    set_vertattr(prog, joints_id, 0, ym::zero4i);
}

//
// This is a public API. See above for documentation.
//
void draw_elem(uint prog, int num, uint bid, etype etype) {
    assert(check_error());
    if (num <= 0) return;
    draw_elems(num, bid, etype);
    assert(check_error());
}

}  // namespace stdshader

}  // namespace yglu
