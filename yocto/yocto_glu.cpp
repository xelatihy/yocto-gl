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

// clang-format off
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl.h>
#include <OpenGL/gl3.h>
#endif
#ifndef YGLU_NO_GLFW
#include <GLFW/glfw3.h>
#endif
// clang-format on

namespace yglu {

//
// Checks for GL error and then prints
//
YGLU_API bool check_error(bool print) {
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
        case GL_STACK_UNDERFLOW: printf("GL_STACK_UNDERFLOW\n"); break;
        case GL_STACK_OVERFLOW: printf("GL_STACK_OVERFLOW\n"); break;
        default: printf("<UNKNOWN GL ERROR>\n"); break;
    }
    return false;
}

//
// Clear window
//
YGLU_API void clear_buffers(const float4& background) {
    glClearColor(background[0], background[1], background[2], background[3]);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

//
// Enable/disable depth test
//
YGLU_API void enable_depth_test(bool enabled) {
    if (enabled)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
}

//
// Enable/disable culling
//
YGLU_API void enable_culling(bool enabled) {
    if (enabled)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
}

//
// Enable/disable wireframe
//
YGLU_API void enable_wireframe(bool enabled) {
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

//
// Enable/disable edges. Attempts to avoid z-fighting but the method is not
// robust.
//
YGLU_API void enable_edges(bool enabled, float tolerance) {
    if (enabled) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDepthRange(0, tolerance);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDepthRange(0, 1);
    }
}

//
// Line width
//
YGLU_API void line_width(float w) { glLineWidth(w); }

//
// Set viewport
//
YGLU_API void set_viewport(const int4& v) {
    glViewport(v[0], v[1], v[2], v[3]);
}

namespace legacy {
//
// This is a public API. See above for documentation.
//
YGLU_API void draw_image(GLuint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom) {
    assert(tid != 0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBindTexture(GL_TEXTURE_2D, tid);
    glEnable(GL_TEXTURE_2D);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, win_w, win_h, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glColor4f(1, 1, 1, 1);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(ox, oy);
    glTexCoord2f(0, 1);
    glVertex2f(ox, oy + zoom * img_h);
    glTexCoord2f(1, 1);
    glVertex2f(ox + zoom * img_w, oy + zoom * img_h);
    glTexCoord2f(1, 0);
    glVertex2f(ox + zoom * img_w, oy);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void read_imagef(float* pixels, int w, int h, int nc) {
    GLuint formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
}

//
// Implementation of make_texture.
//
YGLU_API uint _make_texture(int w, int h, int nc, const void* pixels,
    GLuint type, bool linear, bool mipmap) {
    GLuint formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }
    glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0, formats[nc - 1],
        type, pixels);
    assert(glGetError() == GL_NO_ERROR);
    return id;
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_texture(
    int w, int h, int nc, const float* pixels, bool linear, bool mipmap) {
    return _make_texture(w, h, nc, pixels, GL_FLOAT, linear, mipmap);
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap) {
    return _make_texture(w, h, nc, pixels, GL_UNSIGNED_BYTE, linear, mipmap);
}

//
// Implementation of update_texture.
//
static void _update_texture(uint id, int w, int h, int nc, const void* pixels,
    GLuint type, bool mipmap) {
    GLuint formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(
        GL_TEXTURE_2D, 0, 0, 0, w, h, formats[nc - 1], type, pixels);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void update_texture(
    uint id, int w, int h, int nc, const float* pixels, bool mipmap) {
    _update_texture(id, w, h, nc, pixels, GL_FLOAT, mipmap);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void update_texture(
    uint id, int w, int h, int nc, const unsigned char* pixels, bool mipmap) {
    _update_texture(id, w, h, nc, pixels, GL_UNSIGNED_BYTE, mipmap);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void clear_texture(uint* tid) {
    glDeleteTextures(1, tid);
    *tid = 0;
}

//
// This is a public API. See above for documentation.
//
YGLU_API void begin_frame(const float4x4& camera_xform,
    const float4x4& camera_xform_inv, const float4x4& camera_proj,
    bool eyelight, bool scale_kx) {
    glPushMatrix();
    glPushAttrib(GL_LIGHTING_BIT);

    glEnable(GL_LIGHTING);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
    // glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

#if 0
    if (light_amb != float3{0, 0, 0}) {
        float amb[] = {light_amb[0], light_amb[1], light_amb[2], 1};
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb);
    }
#endif

    if (eyelight) {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        float ke[] = {(scale_kx) ? 3.1415926536f : 1,
            (scale_kx) ? 3.1415926536f : 1, (scale_kx) ? 3.1415926536f : 1, 1};
        float pos[] = {0, 0, 1, 0};
        glLightfv(GL_LIGHT0, GL_DIFFUSE, ke);
        glLightfv(GL_LIGHT0, GL_SPECULAR, ke);
        glLightfv(GL_LIGHT0, GL_POSITION, pos);
    }

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(&camera_proj[0][0]);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(&camera_xform_inv[0][0]);
}

//
// Ends a frame.
//
YGLU_API void end_frame() {
    glPopAttrib();
    glPopMatrix();
}

//
// This is a public API. See above for documentation.
//
YGLU_API void set_lights(const float3& amb, int num, const float3* pos,
    const float3* ke, const ltype* ltype) {
    if (amb != float3{0, 0, 0}) {
        float amb_[] = {amb[0], amb[1], amb[2], 1};
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb_);
    }

    for (auto i = 0; i < std::min(num, 16); i++) {
        glEnable(GL_LIGHT0 + i);
        float ke_[] = {ke[i][0], ke[i][1], ke[i][2], 1};
        float pos_[] = {pos[i][0], pos[i][1], pos[i][2],
            (ltype[i] == ltype::point) ? 1.0f : 0.0f};
        glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, ke_);
        glLightfv(GL_LIGHT0 + i, GL_SPECULAR, ke_);
        glLightfv(GL_LIGHT0 + i, GL_POSITION, pos_);
        if (ltype[i] == ltype::point) {
            glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 0);
            glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0);
            glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 1);
        } else {
            glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 1);
            glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0);
            glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 0);
        }
    }
    for (auto i = num; i < 16; i++) { glDisable(GL_LIGHT0 + i); }

#if 0
    if(light_amb != float3{0,0,0}) {
        glEnable(GL_LIGHT0);
        float amb[] = {light_amb[0], light_amb[1], light_amb[2], 1};
        glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
    }
#endif
}

//
// This is a public API. See above for documentation.
//
YGLU_API void begin_shape(const float4x4& xform) {
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glMultMatrixf((float*)xform.data());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void end_shape() {
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

//
// This is a public API. See above for documentation.
//
YGLU_API float specular_roughness_to_exponent(float rs) {
    return (rs) ? 2 / (rs * rs) - 2 : 1e6;
}

//
// This is a public API. See above for documentation.
//
YGLU_API void set_material(const float3& ke, const float3& kd, const float3& ks,
    float ns, int kd_txt, bool scale_kx) {
    float kds = (scale_kx) ? 1 / 3.1415926536f : 1;
    float kss = (scale_kx) ? (ns + 2) / (2 * 3.1415926536f) : 1;
    // float kss = (scale_kx) ? 1 / 3.1415926536f : 1;
    float ke_[] = {ke[0], ke[1], ke[2], 1};
    float kd_[] = {kd[0] * kds, kd[1] * kds, kd[2] * kds, 1};
    float ks_[] = {ks[0] * kss, ks[1] * kss, ks[2] * kss, 1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, ke_);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, kd_);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ks_);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, ns);
    if (kd_txt >= 0) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, kd_txt);
    } else {
        glDisable(GL_TEXTURE_2D);
    }
}

//
// This is a public API. See above for documentation.
//
YGLU_API float specular_roughness_to_exponent(float r);

//
// This is a public API. See above for documentation.
//
YGLU_API void draw_elems(int num, const int* elem, etype etype,
    const float3* pos, const float3* norm, const float2* texcoord,
    const float3* color) {
    if (!num) return;
    assert(elem);
    auto nc = 0;
    switch (etype) {
        case etype::point:
            glBegin(GL_POINTS);
            nc = 1;
            break;
        case etype::line:
            glBegin(GL_LINES);
            nc = 2;
            break;
        case etype::triangle:
            glBegin(GL_TRIANGLES);
            nc = 3;
            break;
        case etype::quad:
            glBegin(GL_QUADS);
            nc = 4;
            break;
    }
    for (auto j = 0; j < num; j++) {
        auto e = elem + j * nc;
        for (auto i = 0; i < nc; i++) {
            if (norm)
                glNormal3fv((float*)norm + e[i] * 3);
            else
                glNormal3f(0, 0, 1);
            if (texcoord)
                glTexCoord2fv((float*)texcoord + e[i] * 2);
            else
                glTexCoord2f(0, 0);
            if (color)
                glColor3fv((float*)color + e[i] * 3);
            else
                glColor3f(1, 1, 1);
            glVertex3fv((float*)pos + e[i] * 3);
        }
    }
    glEnd();
}
YGLU_API void draw_points(int num, const int* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color) {
    return draw_elems(num, elem, etype::point, pos, norm, texcoord, color);
}
YGLU_API void draw_lines(int num, const int2* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color) {
    return draw_elems(
        num, (const int*)elem, etype::line, pos, norm, texcoord, color);
}
YGLU_API void draw_triangles(int num, const int3* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color) {
    return draw_elems(
        num, (const int*)elem, etype::triangle, pos, norm, texcoord, color);
}

}  // namespace

namespace modern {

//
// This is a public API. See above for documentation.
//
YGLU_API void shade_image(uint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom) {
    shade_image(tid, img_w, img_h, win_w, win_h, ox, oy, zoom, 0, 1);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void shade_image(uint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom, float exposure, float gamma_) {
    static const std::string& vert =
        ""
        "#version 330\n"
        "\n"
        "layout(location = 0) in vec2 vert_texcoord;"
        "uniform vec2 offset;\n"
        "uniform float zoom;\n"
        "uniform vec2 size;\n"
        "uniform vec2 win_size;\n"
        "out vec2 texcoord;"
        "\n"
        "void main() {\n"
        "    texcoord = vert_texcoord.xy;\n"
        "    vec2 pos = offset + size * vert_texcoord.xy * zoom;\n"
        "    vec2 upos = 2 * pos / win_size - vec2(1,1);\n"
        "    upos.y = - upos.y;\n"
        "    gl_Position = vec4(upos.x, upos.y, 0, 1);\n"
        "}\n"
        "";
    static const std::string& frag =
        ""
        "#version 330\n"
        "\n"
        "in vec2 texcoord;\n"
        "out vec4 color;\n"
        "uniform float exposure;\n"
        "uniform float gamma;\n"
        "\n"
        "uniform sampler2D img;\n"
        "\n"
        "void main() {\n"
        "     vec4 c = texture(img,texcoord);\n"
        "     c.xyz = pow(c.xyz*pow(2,exposure),vec3(1/gamma));\n"
        "     color = c;\n"
        "}\n"
        "";

    assert(check_error());

    static uint vao_id = 0;
    if (!vao_id) vao_id = make_vertex_arrays();
    glBindVertexArray(vao_id);

    assert(check_error());

    static uint prog_id = 0;
    if (!prog_id) prog_id = make_program(vert, frag, 0, 0);

    assert(check_error());

    static GLuint tvbo_id = 0;
    if (!tvbo_id) {
        float texcoord[] = {0, 0, 0, 1, 1, 1, 1, 0};
        glGenBuffers(1, &tvbo_id);
        glBindBuffer(GL_ARRAY_BUFFER, tvbo_id);
        glBufferData(
            GL_ARRAY_BUFFER, sizeof(texcoord), texcoord, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    assert(check_error());

    static GLuint evbo_id = 0;
    if (!evbo_id) {
        int elems[] = {0, 1, 2, 0, 2, 3};
        glGenBuffers(1, &evbo_id);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, evbo_id);
        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER, sizeof(elems), elems, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    assert(tid != 0);

    assert(check_error());

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(prog_id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tid);

    float offset[] = {ox, oy};
    float size[] = {(float)img_w, (float)img_h};
    float win_size[] = {(float)win_w, (float)win_h};
    set_uniform(prog_id, "offset", offset, 2, 1);
    set_uniform(prog_id, "zoom", &zoom, 1, 1);
    set_uniform(prog_id, "size", size, 2, 1);
    set_uniform(prog_id, "win_size", win_size, 2, 1);
    set_uniform(prog_id, "exposure", &exposure, 1, 1);
    set_uniform(prog_id, "gamma", &gamma_, 1, 1);
    set_uniform_texture(prog_id, "img", "", tid, 0);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, tvbo_id);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, evbo_id);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glUseProgram(0);
    glDisableVertexAttribArray(0);

    glDisable(GL_BLEND);

    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void read_imagef(float* pixels, int w, int h, int nc) {
    int formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
}

//
// Implementation of make_texture.
//
YGLU_API uint _make_texture(int w, int h, int nc, const void* pixels,
    GLuint type, bool linear, bool mipmap, bool as_float, bool as_srgb) {
    assert(!as_srgb || !as_float);
    int formats_ub[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    int formats_sub[4] = {
        GL_SLUMINANCE, GL_SLUMINANCE_ALPHA, GL_SRGB, GL_SRGB_ALPHA};
    int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
    int* formats =
        (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
    uint id;
    assert(check_error());
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0, formats_ub[nc - 1],
        type, pixels);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    if (mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
            (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }
    assert(check_error());
    return id;
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_texture(int w, int h, int nc, const float* pixels,
    bool linear, bool mipmap, bool as_float) {
    return _make_texture(
        w, h, nc, pixels, GL_FLOAT, linear, mipmap, as_float, false);
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap, bool as_srgb) {
    return _make_texture(
        w, h, nc, pixels, GL_UNSIGNED_BYTE, linear, mipmap, false, as_srgb);
}

//
// Implementation of update_texture.
//
static inline void _update_texture(uint id, int w, int h, int nc,
    const void* pixels, GLuint type, bool mipmap) {
    int formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(
        GL_TEXTURE_2D, 0, 0, 0, w, h, formats[nc - 1], type, pixels);
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void update_texture(
    uint id, int w, int h, int nc, const float* pixels, bool mipmap) {
    _update_texture(id, w, h, nc, pixels, GL_FLOAT, mipmap);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void update_texture(
    uint id, int w, int h, int nc, const unsigned char* pixels, bool mipmap) {
    _update_texture(id, w, h, nc, pixels, GL_UNSIGNED_BYTE, mipmap);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void clear_texture(uint* tid) {
    glDeleteTextures(1, tid);
    *tid = 0;
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_buffer(
    int num, int size, const void* values, bool elements, bool dynamic) {
    auto target = (elements) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    auto bid = (GLuint)0;
    glGenBuffers(1, &bid);
    glBindBuffer(target, bid);
    glBufferData(target, size * num, values,
        (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    glBindBuffer(target, 0);
    return bid;
}

//
// This is a public API. See above for documentation.
//
YGLU_API void update_buffer(uint bid, int num, int size, const void* values,
    bool elements, bool dynamic) {
    auto target = (elements) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    glBindBuffer(target, bid);
    glBufferSubData(target, 0, size * num, values);
    glBindBuffer(target, 0);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void clear_buffer(uint* bid) {
    glDeleteBuffers(1, bid);
    *bid = 0;
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_vertex_arrays() {
    auto aid = (uint)0;
    glGenVertexArrays(1, &aid);
    return aid;
}

//
// This is a public API. See above for documentation.
//
YGLU_API void clear_vertex_arrays(uint* aid) {
    glDeleteVertexArrays(1, aid);
    *aid = 0;
}

//
// This is a public API. See above for documentation.
//
YGLU_API uint make_program(const std::string& vertex,
    const std::string& fragment, uint* vid, uint* fid) {
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
    glBindAttribLocation(gl_program, 0, "vert_pos");
    glBindAttribLocation(gl_program, 1, "vert_norm");
    glBindAttribLocation(gl_program, 2, "vert_texcoord");
    glLinkProgram(gl_program);
    glValidateProgram(gl_program);
    glGetProgramiv(gl_program, GL_LINK_STATUS, &errflags[0]);
    glGetProgramiv(gl_program, GL_VALIDATE_STATUS, &errflags[1]);
    if (!errflags[0] || !errflags[1]) {
        glGetProgramInfoLog(gl_program, 10000, 0, errbuf);
        printf("program not linked\n\n%s\n\n", errbuf);
        return 0;
    }
    assert(glGetError() == GL_NO_ERROR);

    // returns
    if (vid) *vid = gl_shader[0];
    if (fid) *fid = gl_shader[1];
    return gl_program;
}

//
// This is a public API. See above for documentation.
//
YGLU_API void clear_program(uint* pid, uint* vid, uint* fid) {
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
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_uniform(
    uint prog, const std::string& var, const int* val, int nc, int count) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetUniformLocation(prog, var.c_str());
    if (pos < 0) return false;
    switch (nc) {
        case 1: glUniform1iv(pos, count, val); break;
        case 2: glUniform2iv(pos, count, val); break;
        case 3: glUniform3iv(pos, count, val); break;
        case 4: glUniform4iv(pos, count, val); break;
        default: assert(false);
    }
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_uniform(
    uint prog, const std::string& var, const float* val, int nc, int count) {
    assert((nc >= 1 && nc <= 4) || (nc == 16) || (nc == 12));
    int pos = glGetUniformLocation(prog, var.c_str());
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
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_uniform_texture(uint prog, const std::string& var,
    const std::string& varon, uint tid, uint tunit) {
    int pos = glGetUniformLocation(prog, var.c_str());
    int onpos =
        (!varon.empty()) ? glGetUniformLocation(prog, varon.c_str()) : -1;
    if (pos < 0) return false;
    if (tid > 0) {
        glActiveTexture(GL_TEXTURE0 + tunit);
        glBindTexture(GL_TEXTURE_2D, tid);
        glUniform1i(pos, tunit);
        if (onpos >= 0) glUniform1i(onpos, 1);
    } else {
        glActiveTexture(GL_TEXTURE0 + tunit);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(pos, tunit);
        if (onpos >= 0) glUniform1i(onpos, 0);
    }
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_vertattr_ptr(
    uint prog, const std::string& var, const float* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    if (value) {
        glEnableVertexAttribArray(pos);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, value);
    } else {
        glDisableVertexAttribArray(pos);
    }
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_vertattr_buffer(
    uint prog, const std::string& var, uint bid, int nc) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    if (bid) {
        glEnableVertexAttribArray(pos);
        glBindBuffer(GL_VERTEX_ARRAY, bid);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, 0);
        glBindBuffer(GL_VERTEX_ARRAY, 0);
    } else {
        glDisableVertexAttribArray(pos);
    }
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_vertattr_val(
    uint prog, const std::string& var, const float* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var.c_str());
    if (pos < 0) return false;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttrib1fv(pos, value); break;
        case 2: glVertexAttrib2fv(pos, value); break;
        case 3: glVertexAttrib3fv(pos, value); break;
        case 4: glVertexAttrib4fv(pos, value); break;
        default: assert(false); break;
    }
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool set_vertattr(
    uint prog, const std::string& var, uint bid, int nc, const float* def) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var.c_str());
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
    return true;
}

//
// This is a public API. See above for documentation.
//
YGLU_API bool draw_elems(int nelems, uint bid, etype etype) {
    if (!nelems) return true;
    assert(bid);
    assert(check_error());
    int mode = 0, nc = 0;
    switch (etype) {
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
    return true;
}

}  // namespace

namespace stdshader {

//
// This is a public API. See above for documentation.
//
// Filmic tone mapping from
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
//
YGLU_API void make_program(uint* pid, uint* aid) {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    static const std::string& _vert_shader =
        "#version 330\n"
        "\n"
        "layout(location = 0) in vec3 vert_pos;            // vertex position "
        "(in mesh "
        "coordinate frame)\n"
        "layout(location = 1) in vec3 vert_norm;           // vertex normal   "
        "(in mesh "
        "coordinate frame)\n"
        "layout(location = 2) in vec2 vert_texcoord;       // vertex "
        "texcoords\n"
        "layout(location = 3) in vec3 vert_color;          // vertex color\n"
        "\n"
        "uniform mat4 shape_xform;           // shape transform\n"
        "uniform mat4 camera_xform_inv;      // inverse of the camera frame "
        "(as a "
        "matrix)\n"
        "uniform mat4 camera_proj;           // camera projection\n"
        "\n"
        "out vec3 pos;                   // [to fragment shader] vertex "
        "position (in world coordinate)\n"
        "out vec3 norm;                  // [to fragment shader] vertex "
        "normal "
        "(in world coordinate)\n"
        "out vec2 texcoord;              // [to fragment shader] vertex "
        "texture coordinates\n"
        "out vec3 color;                 // [to fragment shader] vertex "
        "color\n"
        "\n"
        "// main function\n"
        "void main() {\n"
        "    pos = (shape_xform * vec4(vert_pos,1)).xyz;\n"
        "    norm = (shape_xform * vec4(vert_norm,0)).xyz;\n"
        "    color = vert_color;\n"
        "    texcoord = vert_texcoord;\n"
        "    gl_Position = camera_proj * camera_xform_inv * shape_xform * "
        "vec4(vert_pos,1);\n"
        "}\n"
        "";

    static const std::string& _frag_shader =
        "#version 330\n"
        "\n"
        "#define TM_DEF 0\n"
        "#define TM_LINEAR 1\n"
        "#define TM_SRGB 2\n"
        "#define TM_GAMMA 3\n"
        "#define TM_FILMIC 4\n"
        "\n"
        "#define pi 3.14159265\n"
        "\n"
        "in vec3 pos;                   // [from vertex shader] position "
        "in "
        "world space\n"
        "in vec3 norm;                  // [from vertex shader] normal in "
        "world space (need normalization)\n"
        "in vec2 texcoord;              // [from vertex shader] texcoord\n"
        "in vec3 color;                 // [from vertex shader] color\n"
        "\n"
        "uniform mat4 camera_xform;          // camera xform\n"
        "\n"
        "uniform vec3 light_amb;             // ambient light\n"
        "uniform int light_num;              // number of lights\n"
        "uniform int light_type[16];         // light type (0 -> point, 1 -> "
        "directional)\n"
        "uniform vec3 light_pos[16];         // light positions\n"
        "uniform vec3 light_ke[16];          // light intensities\n"
        "\n"
        "uniform int material_etype;         // shading surface type\n"
        "uniform vec3 material_ke;           // material ke\n"
        "uniform vec3 material_kd;           // material kd\n"
        "uniform vec3 material_ks;           // material ks\n"
        "uniform float material_rs;          // material rs\n"
        "uniform bool material_use_phong;    // material use phong\n"
        "\n"
        "uniform bool material_txt_ke_on;    // material ke texture on\n"
        "uniform sampler2D material_txt_ke;  // material ke texture\n"
        "uniform bool material_txt_kd_on;    // material kd texture on\n"
        "uniform sampler2D material_txt_kd;  // material kd texture\n"
        "uniform bool material_txt_ks_on;    // material ks texture on\n"
        "uniform sampler2D material_txt_ks;  // material ks texture\n"
        "uniform bool material_txt_rs_on;    // material rs texture on\n"
        "uniform sampler2D material_txt_rs;  // material rs texture\n"
        "\n"
        "uniform float img_exposure;         // image exposure\n"
        "uniform float img_gamma;            // image gamma\n"
        "uniform int img_tonemap;            // image tonemap preset\n"
        "\n"
        "uniform bool shade_eyelight;        // eyelight shading\n"
        "\n"
        "out vec4 frag_color;        // eyelight shading\n"
        "\n"
        "vec3 brdfcos(int et, vec3 kd, vec3 ks, float rs, vec3 n, vec3 wi, "
        "vec3 wo, bool use_phong) {\n"
        "    vec3 wh = normalize(wi+wo);\n"
        "    float ns = 2/(rs*rs)-2;"
        "    float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);"
        "    if(et == 1) {"
        "        return ((1+dot(wo,wi))/2) * kd/pi;\n"
        "    }\n"
        "    if(et == 2) {"
        "        float si = sqrt(1-ndi*ndi), so = sqrt(1-ndo*ndo), sh = "
        "sqrt(1-ndh*ndh);"
        "        if(si <= 0) return vec3(0);"
        "        vec3 diff = si * kd / pi;"
        "        if(sh<=0) return diff;"
        "        float d = ((2+ns)/(2*pi)) * pow(si,ns);"
        "        vec3 spec = si * ks * d / (4*si*so);"
        "        return diff+spec;"
        "    }\n"
        "    if(et == 3 || et == 4) {"
        "        if(ndi<=0 || ndo <=0) return vec3(0);"
        "        vec3 diff = ndi * kd / pi;"
        "        if(ndh<=0) return diff;"
        "        if(use_phong) {"
        "            float d = ((2+ns)/(2*pi)) * pow(ndh,ns);"
        "            vec3 spec = ndi * ks * d / (4*ndi*ndo);"
        "            return diff+spec;"
        "        } else {\n"
        "            float cos2 = ndh * ndh;"
        "            float tan2 = (1 - cos2) / cos2;"
        "            float alpha2 = rs * rs;"
        "            float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * "
        "(alpha2 + tan2));"
        "            float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * "
        "ndo))) / 2;"
        "            float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * "
        "ndi))) / 2;"
        "            float g = 1 / (1 + lambda_o + lambda_i);"
        "            vec3 spec = ndi * ks * d * g / (4*ndi*ndo);"
        "            return diff+spec;"
        "        }\n"
        "    }\n"
        "}\n"
        "\n"
        "vec3 tm_filmic(vec3 x) {"
        "    float a = 2.51f;\n"
        "    float b = 0.03f;\n"
        "    float c = 2.43f;\n"
        "    float d = 0.59f;\n"
        "    float e = 0.14f;\n"
        "    return clamp((x*(a*x+b))/(x*(c*x+d)+e),0,1);\n"
        "}\n"
        "\n"
        "// main\n"
        "void main() {\n"
        "    // view vector\n"
        "    vec3 wo = normalize( (camera_xform*vec4(0,0,0,1)).xyz - pos );"
        "    // re-normalize normals\n"
        "    vec3 n = normalize(norm);\n"
        "    // use faceforward to ensure the normals points toward us\n"
        "    n = faceforward(n,-wo,n);\n"
        "    // get material color from textures\n"
        "    vec3 ke = color * material_ke * ((material_txt_ke_on) ? "
        "texture(material_txt_ke,texcoord).xyz : vec3(1,1,1));\n"
        "    vec3 kd = color * material_kd * ((material_txt_kd_on) ? "
        "texture(material_txt_kd,texcoord).xyz : vec3(1,1,1));\n"
        "    vec3 ks = color * material_ks * ((material_txt_ks_on) ? "
        "texture(material_txt_ks,texcoord).xyz : vec3(1,1,1));\n"
        "    float rs = material_rs * ((material_txt_rs_on) ? "
        "texture(material_txt_rs,texcoord).x : 1);\n"
        "    // emission\n"
        "    vec3 c = ke;\n"
        "    // check early exit\n"
        "    if(kd != vec3(0,0,0) || ks != vec3(0,0,0)) {\n"
        "\n"
        "    if(shade_eyelight) {\n"
        "        vec3 wi = wo;\n"
        "        vec3 wh = normalize(wi+wo);\n"
        "        // accumulate blinn-phong model\n"
        "        c += pi * "
        "brdfcos(material_etype,kd,ks,rs,n,wi,wo,material_use_phong);\n"
        "    } else {\n"
        "        // accumulate ambient\n"
        "        c += light_amb * kd;\n"
        "        // foreach light\n"
        "        for(int i = 0; i < light_num; i ++) {\n"
        "            vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);\n"
        "            if(light_type[i] == 0) {\n"
        "                // compute point light color at pos\n"
        "                cl = light_ke[i] / pow(length(light_pos[i]-pos),2);\n"
        "                // compute light direction at pos\n"
        "                wi = normalize(light_pos[i]-pos);\n"
        "            }\n"
        "            else if(light_type[i] == 1) {\n"
        "                // compute light color\n"
        "                cl = light_ke[i];\n"
        "                // compute light direction\n"
        "                wi = normalize(light_pos[i]);\n"
        "            }\n"
        "            // compute h\n"
        "            vec3 wh = normalize(wi+wo);\n"
        "            // accumulate blinn-phong model\n"
        "            c += cl * "
        "brdfcos(material_etype,kd,ks,rs,n,wi,wo,material_use_phong);\n"
        "        }\n"
        "    }\n"
        "    }\n"
        "\n"
        "    // final color correction\n"
        "    c = c*pow(2,img_exposure);\n"
        "    if(img_tonemap == TM_SRGB || img_tonemap == TM_DEF)"
        "        c = pow(c,vec3(1/2.2));\n"
        "    if(img_tonemap == TM_GAMMA)"
        "        c = pow(c,vec3(1/img_gamma));\n"
        "    if(img_tonemap == TM_FILMIC)"
        "        c = tm_filmic(c);\n"
        "    // output final color by setting gl_FragColor\n"
        "    frag_color = vec4(c,1);\n"
        "}\n"
        "";

    assert(check_error());
    *aid = modern::make_vertex_arrays();
    glBindVertexArray(*aid);
    *pid = modern::make_program(_vert_shader, _frag_shader, 0, 0);
    glBindVertexArray(0);
    assert(check_error());
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif
}

//
// This is a public API. See above for documentation.
//
YGLU_API void begin_frame(uint prog, uint vao, bool shade_eyelight,
    float img_exposure, tonemap_type img_tonemap, float img_gamma,
    const float4x4& camera_xform, const float4x4& camera_xform_inv,
    const float4x4& camera_proj) {
    assert(check_error());
    glUseProgram(prog);
    glBindVertexArray(vao);
    int shade_eyelighti = (shade_eyelight) ? 1 : 0;
    modern::set_uniform(prog, "shade_eyelight", &shade_eyelighti, 1, 1);
    modern::set_uniform(prog, "img_exposure", &img_exposure, 1, 1);
    modern::set_uniform(prog, "img_gamma", &img_gamma, 1, 1);
    auto img_tonemap_int = (int)img_tonemap;
    modern::set_uniform(prog, "img_tonemap", &img_tonemap_int, 1, 1);
    modern::set_uniform(prog, "camera_xform", &camera_xform[0][0], 16, 1);
    modern::set_uniform(
        prog, "camera_xform_inv", &camera_xform_inv[0][0], 16, 1);
    modern::set_uniform(prog, "camera_proj", &camera_proj[0][0], 16, 1);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void end_frame() {
    glBindVertexArray(0);
    glUseProgram(0);
}

//
// This is a public API. See above for documentation.
//
YGLU_API void set_lights(uint prog, const float3& amb, int num, float3* pos,
    float3* ke, ltype* type) {
    assert(check_error());
    modern::set_uniform(prog, "light_amb", &amb[0], 3, 1);
    modern::set_uniform(prog, "light_num", &num, 1, 1);
    modern::set_uniform(prog, "light_pos", (float*)pos, 3, num);
    modern::set_uniform(prog, "light_ke", (float*)ke, 3, num);
    modern::set_uniform(prog, "light_type", (int*)type, 1, num);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void begin_shape(uint prog, const float4x4& xform) {
    assert(check_error());
    modern::set_uniform(prog, "shape_xform", &xform[0][0], 16, 1);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void end_shape() {
    assert(check_error());
    for (int i = 0; i < 16; i++) glDisableVertexAttribArray(i);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void set_material(uint prog, const float3& ke, const float3& kd,
    const float3& ks, float rs, int ke_txt, int kd_txt, int ks_txt, int rs_txt,
    bool use_phong) {
    assert(check_error());
    modern::set_uniform(prog, "material_ke", &ke[0], 3, 1);
    modern::set_uniform(prog, "material_kd", &kd[0], 3, 1);
    modern::set_uniform(prog, "material_ks", &ks[0], 3, 1);
    modern::set_uniform(prog, "material_rs", &rs, 1, 1);
    modern::set_uniform_texture(
        prog, "material_txt_ke", "material_txt_ke_on", ke_txt, 0);
    modern::set_uniform_texture(
        prog, "material_txt_kd", "material_txt_kd_on", kd_txt, 1);
    modern::set_uniform_texture(
        prog, "material_txt_ks", "material_txt_ks_on", ks_txt, 2);
    modern::set_uniform_texture(
        prog, "material_txt_rs", "material_txt_rs_on", rs_txt, 4);
    int use_phongi = use_phong;
    modern::set_uniform(prog, "material_use_phong", &use_phongi, 1, 1);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API float specular_exponent_to_roughness(float n) {
    return std::sqrt(2 / (n + 2));
}

//
// This is a public API. See above for documentation.
//
YGLU_API void set_vert(
    uint prog, uint pos, uint norm, uint texcoord, uint color) {
    assert(check_error());
    float white[3] = {1, 1, 1};
    float zero[3] = {0, 0, 0};
    modern::set_vertattr(prog, "vert_pos", pos, 3, 0);
    modern::set_vertattr(prog, "vert_norm", norm, 3, zero);
    modern::set_vertattr(prog, "vert_texcoord", texcoord, 2, zero);
    modern::set_vertattr(prog, "vert_color", color, 3, white);
    assert(check_error());
}

//
// This is a public API. See above for documentation.
//
YGLU_API void draw_elem(uint prog, int num, uint bid, etype etype) {
    assert(check_error());
    if (num <= 0) return;
    auto etypei = (int)etype;
    modern::set_uniform(prog, "material_etype", &etypei, 1, 1);
    modern::draw_elems(num, bid, etype);
    assert(check_error());
}

YGLU_API void draw_points(uint prog, int num, uint bid) {
    draw_elem(prog, num, bid, etype::point);
}

YGLU_API void draw_lines(uint prog, int num, uint bid) {
    draw_elem(prog, num, bid, etype::line);
}

YGLU_API void draw_triangles(uint prog, int num, uint bid) {
    draw_elem(prog, num, bid, etype::triangle);
}

}  // namespace

}  // namespace

#ifndef YGLU_NO_GLFW

#include <GLFW/glfw3.h>

#ifndef YGLU_NO_NUKLEAR
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
#include "ext/nuklear.h"
#include "ext/nuklear_glfw_gl2.h"
#include "ext/nuklear_glfw_gl3.h"
#endif

namespace yglu {

namespace ui {
//
// Window
//
struct window {
    GLFWwindow* win = nullptr;
    bool legacy_gl = false;
    void* user_pointer = nullptr;
    nk_context* nk_ctx = nullptr;

    int widget_width = 256;

    text_callback text_cb = nullptr;
    refresh_callback refresh_cb = nullptr;
};

//
// Support
//
static inline void _glfw_error_cb(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

//
// Support
//
static inline void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->nk_ctx) {
        nk_glfw3_gl3_char_callback(win->win, key);
        if (nk_item_is_any_active(win->nk_ctx)) return;
    }
    if (win->text_cb) win->text_cb(win, key);
}

//
// Support
//
static inline void _glfw_refresh_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->refresh_cb) win->refresh_cb(win);
}

//
// initialize glfw
//
window* init_window(int width, int height, const std::string& title,
    bool legacy_gl, void* user_pointer) {
    // window
    auto win = new window();
    win->user_pointer = user_pointer;
    win->legacy_gl = legacy_gl;

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

    win->win = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->win);
    glfwSetWindowUserPointer(win->win, win);

    glfwSetErrorCallback(_glfw_error_cb);

// init gl extensions
#ifndef __APPLE__
    if (!glewInit()) return nullptr;
#endif

    return win;
}

//
// initialize glfw
//
void set_callbacks(
    window* win, text_callback text_cb, refresh_callback refresh_cb) {
    win->text_cb = text_cb;
    win->refresh_cb = refresh_cb;
    if (text_cb) glfwSetCharCallback(win->win, _glfw_text_cb);
}

//
// Clear glfw
//
void clear_window(window* win) {
    glfwDestroyWindow(win->win);
    glfwTerminate();
}

//
// Gets the user poiner
void* get_user_pointer(window* win) { return win->user_pointer; }
//

//
// Set window title
//
void set_window_title(window* win, const std::string& title) {
    glfwSetWindowTitle(win->win, title.c_str());
}

//
// Wait events
//
void wait_events(window* win) { glfwWaitEvents(); }

//
// Poll events
//
void poll_events(window* win) { glfwPollEvents(); }

//
// Swap buffers
//
void swap_buffers(window* win) { glfwSwapBuffers(win->win); }

//
// Should close
//
bool should_close(window* win) { return glfwWindowShouldClose(win->win); }

//
// Mouse button
//
int get_mouse_button(window* win) {
    auto mouse1 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
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
int2 get_mouse_pos(window* win) {
    double x, y;
    glfwGetCursorPos(win->win, &x, &y);
    return {(int)x, (int)y};
}

//
// Mouse position
//
float2 get_mouse_posf(window* win) {
    double x, y;
    glfwGetCursorPos(win->win, &x, &y);
    return {(float)x, (float)y};
}

//
// Window size
//
int2 get_window_size(window* win) {
    auto ret = int2{0, 0};
    glfwGetWindowSize(win->win, &ret[0], &ret[1]);
    return ret;
}

//
// Framebuffer size
//
int2 get_framebuffer_size(window* win) {
    auto ret = int2{0, 0};
    glfwGetFramebufferSize(win->win, &ret[0], &ret[1]);
    return ret;
}

//
// Read pixels
//
std::vector<byte4> get_screenshot(
    window* win, int2& wh, bool flipy, bool back) {
    wh = get_framebuffer_size(win);
    auto pixels = std::vector<byte4>(wh[0] * wh[1]);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(0, 0, wh[0], wh[1], GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    if (flipy) {
        std::vector<byte4> line(wh[0]);
        for (int j = 0; j < wh[1] / 2; j++) {
            memcpy(line.data(), pixels.data() + j * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + j * wh[0] * 4,
                pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, line.data(),
                wh[0] * 4);
        }
    }
    return pixels;
}

#ifndef YGLU_NO_NUKLEAR

//
// Nuklear
//
void init_widgets(window* win) {
    glfwSetScrollCallback(win->win, nk_gflw3_scroll_callback);
    if (win->legacy_gl) {
        win->nk_ctx = nk_glfw3_gl2_init(win->win, NK_GLFW3_GL2_DEFAULT);
        nk_font_atlas* atlas;
        nk_glfw3_gl2_font_stash_begin(&atlas);
        nk_glfw3_gl2_font_stash_end();
    } else {
        win->nk_ctx = nk_glfw3_gl3_init(win->win, NK_GLFW3_GL3_DEFAULT);
        nk_font_atlas* atlas;
        nk_glfw3_gl3_font_stash_begin(&atlas);
        nk_glfw3_gl3_font_stash_end();
    }
}

//
// Nuklear
//
void clear_widgets(window* win) {
    if (win->legacy_gl) {
        nk_glfw3_gl2_shutdown();
    } else {
        nk_glfw3_gl3_shutdown();
    }
}

//
// Begin draw widget
//
bool begin_widgets(window* win) {
    if (!win->nk_ctx) return false;
    auto window_size = get_window_size(win);
    if (win->legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    auto visible = nk_begin(win->nk_ctx, "yshade",
        nk_rect(window_size[0] - win->widget_width, 0, win->widget_width,
            window_size[1]),
        NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
            NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE);
    // if(visible) dynamic_widget_layout(win, 1);
    return visible;
}

//
// Dynamic layout for next widgets
//
void dynamic_widget_layout(window* win, int n) {
    nk_layout_row_dynamic(win->nk_ctx, 30, n);
}

//
// End draw widget
//
void end_widgets(window* win) {
    nk_end(win->nk_ctx);

    if (win->legacy_gl) {
        nk_glfw3_gl2_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    } else {
        nk_glfw3_gl3_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    }
}

//
// Label widget
//
void label_widget(window* win, const std::string& lbl) {
    nk_label(win->nk_ctx, lbl.c_str(), NK_TEXT_LEFT);
}

//
// Label and int widget
//
void int_label_widget(window* win, const std::string& lbl, int val) {
    nk_value_int(win->nk_ctx, lbl.c_str(), val);
}

//
// Label and float widget
//
void float_label_widget(window* win, const std::string& lbl, float val) {
    nk_value_float(win->nk_ctx, lbl.c_str(), val);
}

//
// Label widget
//
void int_widget(
    window* win, const std::string& lbl, int* val, int min, int max, int incr) {
    nk_property_int(win->nk_ctx, lbl.c_str(), min, val, max, incr, incr);
}

//
// Label widget
//
void float_widget(window* win, const std::string& lbl, float* val, float min,
    float max, float incr) {
    nk_property_float(win->nk_ctx, lbl.c_str(), min, val, max, incr, incr);
}

//
// Enum widget
//
void enum_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels) {
    const char* items[100];
    auto pos = -1;
    for (auto i = 0; i < labels.size(); i++) {
        items[i] = labels[i].first.c_str();
        if (labels[i].second == *val) pos = i;
    }
    nk_combobox(win->nk_ctx, items, (int)labels.size(), &pos, 20, {100, 100});
    *val = labels[pos].second;
}

//
// Bool widget
//
void bool_widget(window* win, const std::string& lbl, bool* val) {
    *val = nk_check_label(win->nk_ctx, lbl.c_str(), *val);
}

//
// Button widget
//
bool button_widget(window* win, const std::string& lbl) {
    return nk_button_label(win->nk_ctx, lbl.c_str());
}

//
// Whether widget are active
//
bool get_widget_active(window* win) {
    if (!win->nk_ctx) return false;
    return nk_item_is_any_active(win->nk_ctx);
}
#endif

}  // namespace

}  // namespace

#endif
