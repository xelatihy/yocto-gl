//
// YOCTO_GLU: a set of utilities to draw on screen with OpenGL 2.1. Mostly
// used to make quick viewers. Not sure it is helpful to others (OpenGL is
// really a portability mess, but there is nothing else).
//

//
// FEATURES:
//
// 0. include this file (more compilation options below)
// 1. image viewing
// - yg_draw_image(textureid, image size, window size, offset, zoom)
// - with exposure/gamma yg_shade_image(as above, exposure gamma)
// 2. texture utilies to quickly create/update textures
// - textureid = yg_make_texture(image data, load as floats, filter on/off)
// - yg_update_texture(textureid, image data)
// - yg_clear_texture(textureid)
// 3. program utilities
// - yg_make_program(vertex code, fragment code, output ids)
// - yg_clear_program(ids)
// - yg_set_uniformi(progid, varname, int values, counts)
// - yg_set_uniformf(progid, varname, float values, counts)
// - yg_set_uniformt(progid, varname, textureid, textureunit)
// - yg_set_vertattr(progid, varname, values, count)
// - yg_set_vertattr_ptr(progid, varname, values, count)
// - yg_set_vertattr_val(progid, varname, values, count)
// - yg_draw_elems(element data)
// 4. a standard shader for GGX fragment shading and multiple lights
// - yg_stdshader_make_program()
// - yg_stdshader_begin_frame(image and camera params)
// - yg_stdshader_set_lights(light params)
// - yg_stdshader_begin_shape(shape xform)
// - yg_stdshader_end_shape()
// - yg_stdshader_set_material(material and texture data)
// - yg_stdshader_set_vert(vertex data)
// - yg_stdshader_draw_elem(element data)
// - yg_stdshader_end_frame()
//
// The interface for each function is described in details in the interface
// section of this file.

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use in
// C/C++.
// To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.1: C/C++ implementation
// - v 0.0: initial release in C99
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YG_H_
#define _YG_H_

// compilation options
#ifdef __cplusplus
#ifndef YGL_DECLARATION
#define YGL_API inline
#define YGLC_API inline
#else
#define YGL_API
#define YGLC_API extern "C"
#endif
#endif

#ifndef __cplusplus
#define YGLC_API extern
#include <stdbool.h>
#endif

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Shape types
//
enum {
    yg_etype_point = 1,     // points
    yg_etype_line = 2,      // lines
    yg_etype_triangle = 3,  // triangles
    yg_etype_quad = 4       // quads
};

//
// Light types
//
enum {
    yg_ltype_point = 0,       // point lights
    yg_ltype_directional = 1  // directional lights
};

// IMAGE FUNCTIONS -------------------------------------------------------------

//
// Draw an texture tid of size img_w, img_h on a window of size win_w, win_h
// with top-left corner at ox, oy with a zoom zoom.
//
YGLC_API void yg_draw_image(int tid, int img_w, int img_h, int win_w, int win_h,
                            float ox, float oy, float zoom);

//
// As above but includes an exposure/gamma correction.
//
YGLC_API void yg_shade_image(int tid, int img_w, int img_h, int win_w,
                             int win_h, float ox, float oy, float zoom,
                             float exposure, float gamma_);

//
// Reads an image back to memory.
//
YGLC_API void yg_read_imagef(float* pixels, int w, int h, int nc);

// TEXTURE FUNCTIONS -----------------------------------------------------------

//
// Creates a texture with pixels values of size w, h with nc number of
// components (1-4). Internall use float if as_float and filtering if filter.
// Returns the texture id.
//
YGLC_API int yg_make_texture(const float* pixels, int w, int h, int nc,
                             bool as_float, bool filter);

//
// Updates the texture tid with new image data.
//
YGLC_API int yg_update_texture(int tid, const float* pixels, int w, int h,
                               int nc);

//
// Destroys the texture tid.
//
YGLC_API void yg_clear_texture(int* tid);

// PROGRAM FUNCTIONS -----------------------------------------------------------

//
// Creates and OpenGL program from vertex and fragment code. Returns the
// program id. Optionally return vertex and fragment shader ids.
//
YGLC_API int yg_make_program(const char* vertex, const char* fragment, int* vid,
                             int* fid);

//
// Destroys the program pid and optionally the sahders vid and fid.
//
YGLC_API void yg_clear_program(int* pid, int* vid, int* fid);

//
// Set uniform integer values val for program prog and variable var.
// The values have nc number of components (1-4) and count elements
// (for arrays).
//
YGLC_API int yg_set_uniformi(int prog, const char* var, const int* val, int nc,
                             int count);

//
// Set uniform float values val for program prog and variable var.
// The values have nc number of components (1-4) and count elements
// (for arrays).
//
YGLC_API int yg_set_uniformf(int prog, const char* var, const float* val,
                             int nc, int count);

//
// Set uniform texture id tid and unit tunit for program prog and variable var.
// Optionally sets the int variable varon to 0/1 whether the texture is enable
// on not.
//
YGLC_API int yg_set_uniformt(int prog, const char* var, const char* varon,
                             int tid, int tunit);

//
// Sets a vartex attribute pointer for program prog and variable var.
// The attribute has nc components and per-vertex values values.
//
YGLC_API int yg_set_vertattr_ptr(int prog, const char* var, const float* value,
                                 int nc);

//
// Sets a constant value for a vertex attribute for program prog and
// variable var. The attribute has nc components.
//
YGLC_API int yg_set_vertattr_val(int prog, const char* var, const float* value,
                                 int nc);

//
// Sets a vartex attribute for program prog and variable var. The attribute
// has nc components and either per-vertex values value or a single value def
// (if value is NULL). Convenience wraper to above functions.
//
YGLC_API int yg_set_vertattr(int prog, const char* var, const float* value,
                             int nc, const float* def);

//
// Draws nelems elements elem of type etype.
//
YGLC_API int yg_draw_elems(int nelems, const int* elem, int etype);

// STANDARD SHADER FUNCTIONS ---------------------------------------------------

//
// Initialize a standard shader.
//
YGLC_API int yg_stdshader_make_program();

//
// Starts a frame by setting exposure/gamma values, camera transforms and
// projection. Sets also whether to use full shading or a quick eyelight
// preview.
//
YGLC_API void yg_stdshader_begin_frame(int prog, bool shade_eyelight,
                                       float img_exposure, float img_gamma,
                                       const float camera_xform[16],
                                       const float camera_xform_inv[16],
                                       const float camera_proj[16]);

//
// Ends a frame.
//
YGLC_API void yg_stdshader_end_frame();

//
// Set num lights with position pos, color ke, type ltype. Also set the ambient
// illumination amb.
//
YGLC_API void yg_stdshader_set_lights(int prog, const float amb[3], int num,
                                      float* pos, float* ke, int* ltype);

//
// Begins drawing a shape with transform xform.
//
YGLC_API void yg_stdshader_begin_shape(int prog, const float xform[16]);

//
// End shade drawing.
//
YGLC_API void yg_stdshader_end_shape();

//
// Set material values with emission ke, diffuse kd, specular ks and
// specular roughness rs. Indicates textures ids with the correspoinding XXX_txt
// variables. Uses GGX by default, but can switch to Phong is needed. Works
// for points, lines and triangle/quads based on etype.
//
YGLC_API void yg_stdshader_set_material(int prog, int etype, float ke[3],
                                        float kd[3], float ks[3], float rs,
                                        int ke_txt, int kd_txt, int ks_txt,
                                        int rs_txt, bool use_phong);

//
// Convertes a phong exponent to roughness.
//
YGLC_API float yg_specular_exponent_to_roughness(float n);

//
// Set vertex data with position pos, normals norm, texture coordinates texcoord
// and per-vertex color color.
//
YGLC_API void yg_stdshader_set_vert(int prog, float* pos, float* norm,
                                    float* texcoord, float* color);

//
// Draw num elements elem of type etype.
//
YGLC_API void yg_stdshader_draw_elem(int prog, int num, int* elem, int etype);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION)

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#endif

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_draw_image(int tid, int img_w, int img_h, int win_w, int win_h,
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
YGLC_API void yg_shade_image(int tid, int img_w, int img_h, int win_w,
                             int win_h, float ox, float oy, float zoom,
                             float exposure, float gamma_) {
    static const char* vert =
        ""
        "#version 120\n"
        "\n"
        "varying vec2 texcoord;"
        "\n"
        "void main() {\n"
        "    texcoord = gl_MultiTexCoord0.xy;\n"
        "    gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;\n"
        "}\n"
        "";
    static const char* frag =
        ""
        "#version 120\n"
        "\n"
        "varying vec2 texcoord;\n"
        "uniform float exposure;\n"
        "uniform float gamma;\n"
        "\n"
        "uniform sampler2D img;\n"
        "\n"
        "void main() {\n"
        "     vec4 c = texture2D(img,texcoord);\n"
        "     c.xyz = pow(c.xyz*pow(2,exposure),vec3(1/gamma));\n"
        "     gl_FragColor = c;\n"
        "}\n"
        "";

    static int prog_id = 0;
    if (!prog_id) {
        prog_id = yg_make_program(vert, frag, 0, 0);
    }

    assert(tid != 0);

    glUseProgram(prog_id);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tid);
    yg_set_uniformf(prog_id, "exposure", &exposure, 1, 1);
    yg_set_uniformf(prog_id, "gamma", &gamma_, 1, 1);
    yg_set_uniformt(prog_id, "img", 0, tid, 0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, win_w, win_h, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
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
    glDisable(GL_BLEND);
    glUseProgram(0);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_read_imagef(float* pixels, int w, int h, int nc) {
    int formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_make_texture(const float* pixels, int w, int h, int nc,
                             bool as_float, bool filter) {
    int formats_ub[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    int formats_f[4] = {GL_LUMINANCE32F_ARB, GL_LUMINANCE_ALPHA32F_ARB,
                        GL_RGB32F_ARB, GL_RGBA32F_ARB};
    int* formats = (as_float) ? formats_f : formats_ub;
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    if (filter) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    } else {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }
    glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0, formats_ub[nc - 1],
                 GL_FLOAT, pixels);
    assert(glGetError() == GL_NO_ERROR);
    return id;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_update_texture(int id, const float* pixels, int w, int h,
                               int nc) {
    int formats[4] = {GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, formats[nc - 1], GL_FLOAT,
                    pixels);
    return id;
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_clear_texture(int* tid) {
    glDeleteTextures(1, (GLuint*)tid);
    *tid = 0;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_make_program(const char* vertex, const char* fragment, int* vid,
                             int* fid) {
    int errflags[2];
    char errbuf[10000];
    int gl_shader[2] = {0, 0};
    const char* code[2] = {vertex, fragment};
    int type[2] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER};
    for (int i = 0; i < 2; i++) {
        gl_shader[i] = glCreateShader(type[i]);
        glShaderSource(gl_shader[i], 1, (const char**)&code[i], NULL);
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
YGLC_API void yg_clear_program(int* pid, int* vid, int* fid) {
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
YGLC_API int yg_set_uniformi(int prog, const char* var, const int* val, int nc,
                             int count) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetUniformLocation(prog, var);
    if (pos < 0) return 0;
    switch (nc) {
        case 1: glUniform1iv(pos, count, val); break;
        case 2: glUniform2iv(pos, count, val); break;
        case 3: glUniform3iv(pos, count, val); break;
        case 4: glUniform4iv(pos, count, val); break;
        default: assert(false);
    }
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_set_uniformf(int prog, const char* var, const float* val,
                             int nc, int count) {
    assert((nc >= 1 && nc <= 4) || (nc == 16) || (nc == 12));
    int pos = glGetUniformLocation(prog, var);
    if (pos < 0) return 0;
    switch (nc) {
        case 1: glUniform1fv(pos, count, val); break;
        case 2: glUniform2fv(pos, count, val); break;
        case 3: glUniform3fv(pos, count, val); break;
        case 4: glUniform4fv(pos, count, val); break;
        case 12: glUniformMatrix4x3fv(pos, count, false, val); break;
        case 16: glUniformMatrix4fv(pos, count, false, val); break;
        default: assert(false); return 0;
    }
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_set_uniformt(int prog, const char* var, const char* varon,
                             int tid, int tunit) {
    int pos = glGetUniformLocation(prog, var);
    int onpos = (varon) ? glGetUniformLocation(prog, varon) : -1;
    if (pos < 0) return 0;
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
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_set_vertattr_ptr(int prog, const char* var, const float* value,
                                 int nc) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var);
    if (pos < 0) return 0;
    if (value) {
        glEnableVertexAttribArray(pos);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, value);
    } else {
        glDisableVertexAttribArray(pos);
    }
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_set_vertattr_val(int prog, const char* var, const float* value,
                                 int nc) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var);
    if (pos < 0) return 0;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttrib1fv(pos, value); break;
        case 2: glVertexAttrib2fv(pos, value); break;
        case 3: glVertexAttrib3fv(pos, value); break;
        case 4: glVertexAttrib4fv(pos, value); break;
        default: assert(false); break;
    }
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_set_vertattr(int prog, const char* var, const float* value,
                             int nc, const float* def) {
    assert(nc >= 1 && nc <= 4);
    int pos = glGetAttribLocation(prog, var);
    if (pos < 0) return 0;
    if (value) {
        glEnableVertexAttribArray(pos);
        glVertexAttribPointer(pos, nc, GL_FLOAT, false, 0, value);
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
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_draw_elems(int nelems, const int* elem, int etype) {
    int mode = 0;
    switch (etype) {
        case yg_etype_point: mode = GL_POINTS; break;
        case yg_etype_line: mode = GL_LINES; break;
        case yg_etype_triangle: mode = GL_TRIANGLES; break;
        case yg_etype_quad: mode = GL_QUADS; break;
        default: assert(false);
    };
    glDrawElements(mode, nelems * etype, GL_UNSIGNED_INT, elem);
    return 1;
}

//
// This is a public API. See above for documentation.
//
YGLC_API int yg_stdshader_make_program() {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    static const char* yg__vert_shader =
        "#version 120\n"
        "\n"
        "attribute vec3 vert_pos;            // vertex position (in mesh "
        "coordinate frame)\n"
        "attribute vec3 vert_norm;           // vertex normal   (in mesh "
        "coordinate frame)\n"
        "attribute vec2 vert_texcoord;       // vertex texcoords\n"
        "attribute vec3 vert_color;          // vertex color\n"
        "\n"
        "uniform mat4 shape_xform;           // shape transform\n"
        "uniform mat4 camera_xform_inv;      // inverse of the camera frame "
        "(as a "
        "matrix)\n"
        "uniform mat4 camera_proj;           // camera projection\n"
        "\n"
        "varying vec3 pos;                   // [to fragment shader] vertex "
        "position (in world coordinate)\n"
        "varying vec3 norm;                  // [to fragment shader] vertex "
        "normal "
        "(in world coordinate)\n"
        "varying vec2 texcoord;              // [to fragment shader] vertex "
        "texture coordinates\n"
        "varying vec3 color;                 // [to fragment shader] vertex "
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

    static const char* yg__frag_shader =
        "#version 120\n"
        "\n"
        "#define pi 3.14159265\n"
        "\n"
        "varying vec3 pos;                   // [from vertex shader] position "
        "in "
        "world space\n"
        "varying vec3 norm;                  // [from vertex shader] normal in "
        "world space (need normalization)\n"
        "varying vec2 texcoord;              // [from vertex shader] texcoord\n"
        "varying vec3 color;                 // [from vertex shader] color\n"
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
        "\n"
        "uniform bool shade_eyelight;        // eyelight shading\n"
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
        "texture2D(material_txt_ke,texcoord).xyz : vec3(1,1,1));\n"
        "    vec3 kd = color * material_kd * ((material_txt_kd_on) ? "
        "texture2D(material_txt_kd,texcoord).xyz : vec3(1,1,1));\n"
        "    vec3 ks = color * material_ks * ((material_txt_ks_on) ? "
        "texture2D(material_txt_ks,texcoord).xyz : vec3(1,1,1));\n"
        "    float rs = material_rs * ((material_txt_rs_on) ? "
        "texture2D(material_txt_rs,texcoord).x : 1);\n"
        "    // emission\n"
        "    vec3 c = ke;\n"
        "    // check early exit\n"
        "    if(kd == vec3(0,0,0) && ks == vec3(0,0,0)) {\n"
        "        c = pow(c*pow(2,img_exposure),vec3(1/img_gamma));\n"
        "        gl_FragColor = vec4(c,1);\n"
        "        return;\n"
        "    }\n"
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
        "\n"
        "    //final color correction\n"
        "    c = pow(c*pow(2,img_exposure),vec3(1/img_gamma));\n"
        "    // output final color by setting gl_FragColor\n"
        "    gl_FragColor = vec4(c,1);\n"
        "}\n"
        "";

    return yg_make_program(yg__vert_shader, yg__frag_shader, 0, 0);
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_begin_frame(int prog, bool shade_eyelight,
                                       float img_exposure, float img_gamma,
                                       const float camera_xform[16],
                                       const float camera_xform_inv[16],
                                       const float camera_proj[16]) {
    glUseProgram(prog);
    int shade_eyelighti = (shade_eyelight) ? 1 : 0;
    yg_set_uniformi(prog, "shade_eyelight", &shade_eyelighti, 1, 1);
    yg_set_uniformf(prog, "img_exposure", &img_exposure, 1, 1);
    yg_set_uniformf(prog, "img_gamma", &img_gamma, 1, 1);
    yg_set_uniformf(prog, "camera_xform", camera_xform, 16, 1);
    yg_set_uniformf(prog, "camera_xform_inv", camera_xform_inv, 16, 1);
    yg_set_uniformf(prog, "camera_proj", camera_proj, 16, 1);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_end_frame() { glUseProgram(0); }

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_set_lights(int prog, const float amb[3], int num,
                                      float* pos, float* ke, int* type) {
    yg_set_uniformf(prog, "light_amb", amb, 3, 1);
    yg_set_uniformi(prog, "light_num", &num, 1, 1);
    yg_set_uniformf(prog, "light_pos", pos, 3, num);
    yg_set_uniformf(prog, "light_ke", ke, 3, num);
    yg_set_uniformi(prog, "light_type", type, 1, num);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_begin_shape(int prog, const float* xform) {
    yg_set_uniformf(prog, "shape_xform", xform, 16, 1);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_end_shape() {
    for (int i = 0; i < 16; i++) glEnableVertexAttribArray(i);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_set_material(int prog, int etype, float ke[3],
                                        float kd[3], float ks[3], float rs,
                                        int ke_txt, int kd_txt, int ks_txt,
                                        int rs_txt, bool use_phong) {
    yg_set_uniformi(prog, "material_etype", &etype, 1, 1);
    yg_set_uniformf(prog, "material_ke", ke, 3, 1);
    yg_set_uniformf(prog, "material_kd", kd, 3, 1);
    yg_set_uniformf(prog, "material_ks", ks, 3, 1);
    yg_set_uniformf(prog, "material_rs", &rs, 1, 1);
    yg_set_uniformt(prog, "material_txt_ke", "material_txt_ke_on", ke_txt, 0);
    yg_set_uniformt(prog, "material_txt_kd", "material_txt_kd_on", kd_txt, 1);
    yg_set_uniformt(prog, "material_txt_ks", "material_txt_ks_on", ks_txt, 2);
    yg_set_uniformt(prog, "material_txt_rs", "material_txt_rs_on", rs_txt, 4);
    int use_phongi = use_phong;
    yg_set_uniformi(prog, "material_use_phong", &use_phongi, 1, 1);
}

//
// This is a public API. See above for documentation.
//
YGLC_API float yg_specular_exponent_to_roughness(float n) {
    return sqrtf(2 / (n + 2));
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_set_vert(int prog, float* pos, float* norm,
                                    float* texcoord, float* color) {
    float white[3] = {1, 1, 1};
    float zero[3] = {0, 0, 0};
    yg_set_vertattr(prog, "vert_pos", pos, 3, 0);
    yg_set_vertattr(prog, "vert_norm", norm, 3, zero);
    yg_set_vertattr(prog, "vert_texcoord", texcoord, 2, zero);
    yg_set_vertattr(prog, "vert_color", color, 3, white);
}

//
// This is a public API. See above for documentation.
//
YGLC_API void yg_stdshader_draw_elem(int prog, int num, int* elem, int etype) {
    yg_draw_elems(num, elem, etype);
}

#endif

#endif
