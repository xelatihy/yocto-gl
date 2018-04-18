//
// Implementation for Yocto/GLUtils.
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

#include "yocto_glutils.h"
#include "yocto_utils.h"

#include <cassert>
#include <unordered_map>

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#define GLFW_INCLUDE_GLCOREARB
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw_gl3.h"
unsigned int imgui_extrafont_compressed_size();
const unsigned int* imgui_extrafont_compressed_data();

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPENGL
// -----------------------------------------------------------------------------
namespace ygl {
// Checks for GL error and then prints
bool check_glerror(bool print) {
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

// Clear window
void clear_glbuffers(const vec4f& background) {
    assert(check_glerror());
    glClearColor(background.x, background.y, background.z, background.w);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    assert(check_glerror());
}

// Enable/disable depth test
void enable_gldepth_test(bool enabled) {
    assert(check_glerror());
    if (enabled) {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
    } else {
        glDisable(GL_DEPTH_TEST);
    }
    assert(check_glerror());
}

// Enable/disable culling
void enable_glculling(bool enabled, bool front, bool back) {
    assert(check_glerror());
    if (enabled && (front || back)) {
        glEnable(GL_CULL_FACE);
        if (front && back)
            glCullFace(GL_FRONT_AND_BACK);
        else if (front)
            glCullFace(GL_FRONT);
        else if (back)
            glCullFace(GL_BACK);
    } else {
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    assert(check_glerror());
}

// Enable/disable wireframe
void enable_glwireframe(bool enabled) {
    assert(check_glerror());
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    assert(check_glerror());
}

// Enable/disable blending
void enable_glblending(bool enabled) {
    assert(check_glerror());
    if (enabled) {
        glEnable(GL_BLEND);
    } else {
        glDisable(GL_BLEND);
    }
    assert(check_glerror());
}

// Set blending to over operator
void set_glblend_over() {
    assert(check_glerror());
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    assert(check_glerror());
}

// Line width
void gl_line_width(float w) {
    assert(check_glerror());
    glLineWidth(min(max(w, 0.0f), 1.0f));
    assert(check_glerror());
}

// Set viewport
void set_glviewport(const vec4i& v) {
    assert(check_glerror());
    glViewport(v.x, v.y, v.z, v.w);
    assert(check_glerror());
}

// Set viewport
void set_glviewport(const vec2i& v) {
    assert(check_glerror());
    glViewport(0, 0, v.x, v.y);
    assert(check_glerror());
}

// This is a public API. See above for documentation.
void read_glimagef(float* pixels, int w, int h, int nc) {
    assert(check_glerror());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of update_texture.
void update_gltexture(gltexture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb) {
    auto refresh = !txt.tid || txt.width != w || txt.height != h ||
                   txt.ncomp != nc || txt.as_float != as_float ||
                   txt.as_srgb != as_srgb || txt.mipmap != mipmap ||
                   txt.linear != linear;
    txt.width = w;
    txt.height = h;
    txt.ncomp = nc;
    txt.as_float = as_float;
    txt.as_srgb = as_srgb;
    txt.mipmap = mipmap;
    txt.linear = linear;
    assert(!as_srgb || !as_float);
    assert(check_glerror());
    if (w * h) {
        int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
        int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
        int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
        int* formats =
            (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
        assert(check_glerror());
        if (!txt.tid) glGenTextures(1, &txt.tid);
        glBindTexture(GL_TEXTURE_2D, txt.tid);
        if (refresh) {
            glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0,
                formats_ub[nc - 1], (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE,
                pixels);
        } else {
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, formats_ub[nc - 1],
                (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE, pixels);
        }
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
    } else {
        if (txt.tid) {
            glBindTexture(GL_TEXTURE_2D, txt.tid);
            glDeleteTextures(1, &txt.tid);
            txt.tid = 0;
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    }
    assert(check_glerror());
}

// Binds a texture to a texture unit
void bind_gltexture(const gltexture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, txt.tid);
}

// Unbinds
void unbind_gltexture(const gltexture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Destroys the texture tid.
void clear_gltexture(gltexture& txt) {
    assert(check_glerror());
    glDeleteTextures(1, &txt.tid);
    txt.tid = 0;
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL VERTEX/ELEMENTS BUFFERS
// -----------------------------------------------------------------------------
namespace ygl {

// Updates the buffer with new data.
void update_glbuffer(glvertex_buffer& buf, bool elems, int num, int ncomp,
    const void* values, bool dynamic) {
    auto resize =
        !buf.bid || num * ncomp != buf.num * buf.ncomp || elems != buf.elems;
    buf.num = num;
    buf.ncomp = ncomp;
    buf.elems = elems;
    auto target = (elems) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    auto esize = (elems) ? sizeof(int) : sizeof(float);
    auto dmode = (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW;
    assert(check_glerror());
    if (num) {
        if (!buf.bid) glGenBuffers(1, &buf.bid);
        glBindBuffer(target, buf.bid);
        if (resize) {
            glBufferData(target, buf.num * buf.ncomp * esize, values, dmode);
        } else {
            glBufferSubData(target, 0, buf.num * buf.ncomp * esize, values);
        }
        glBindBuffer(target, 0);
    } else {
        if (buf.bid) {
            glBindBuffer(target, buf.bid);
            glDeleteBuffers(1, &buf.bid);
            buf.bid = 0;
            glBindBuffer(target, 0);
        }
    }
    assert(check_glerror());
}

// Bind the buffer at a particular attribute location
void bind_glbuffer(const glvertex_buffer& buf, uint vattr) {
    if (buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    glEnableVertexAttribArray(vattr);
    glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
    glVertexAttribPointer(vattr, buf.ncomp, GL_FLOAT, false, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Unbind the buffer
void unbind_glbuffer(const glvertex_buffer& buf, uint vattr) {
    if (buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    glDisableVertexAttribArray(vattr);
}

// Draws elements.
void draw_glelems(const glvertex_buffer& buf) {
    if (!buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    if (!buf.bid) return;
    assert(check_glerror());
    int mode = 0;
    switch (buf.ncomp) {
        case 1: mode = GL_POINTS; break;
        case 2: mode = GL_LINES; break;
        case 3: mode = GL_TRIANGLES; break;
        case 4: mode = GL_QUADS; break;
        default: assert(false);
    };
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glDrawElements(mode, buf.ncomp * buf.num, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(check_glerror());
}

// Unbind the buffer
void unbind_glbuffer(uint vattr) { glDisableVertexAttribArray(vattr); }

// Destroys the buffer
void clear_glbuffer(glvertex_buffer& buf) {
    assert(check_glerror());
    glDeleteBuffers(1, &buf.bid);
    buf.bid = 0;
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Creates and OpenGL program from vertex and fragment code. Returns the
// program id. Optionally return vertex and fragment shader ids. A VAO
// is created.
glprogram make_glprogram(
    const std::string& vertex, const std::string& fragment) {
    auto prog = glprogram();

    assert(check_glerror());
    glGenVertexArrays(1, &prog.vao);
    glBindVertexArray(prog.vao);
    assert(check_glerror());

    int errflags;
    char errbuf[10000];

    // create vertex
    prog.vid = glCreateShader(GL_VERTEX_SHADER);
    const char* vertex_str = vertex.c_str();
    glShaderSource(prog.vid, 1, &vertex_str, NULL);
    glCompileShader(prog.vid);
    glGetShaderiv(prog.vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create fragment
    prog.fid = glCreateShader(GL_FRAGMENT_SHADER);
    const char* fragment_str = fragment.c_str();
    glShaderSource(prog.fid, 1, &fragment_str, NULL);
    glCompileShader(prog.fid);
    glGetShaderiv(prog.fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create program
    prog.pid = glCreateProgram();
    glAttachShader(prog.pid, prog.vid);
    glAttachShader(prog.pid, prog.fid);
    glLinkProgram(prog.pid);
    glValidateProgram(prog.pid);
    glGetProgramiv(prog.pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(prog.pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    assert(check_glerror());

    glBindVertexArray(0);
    assert(check_glerror());

    return prog;
}

// Destroys the program pid and optionally the sahders vid and fid.
void clear_program(glprogram& prog) {
    assert(check_glerror());
    glDetachShader(prog.pid, prog.vid);
    glDeleteShader(prog.vid);
    prog.vid = 0;
    glDetachShader(prog.pid, prog.fid);
    glDeleteShader(prog.fid);
    prog.fid = 0;
    glDeleteProgram(prog.pid);
    prog.pid = 0;
    glDeleteVertexArrays(1, &prog.vao);
    prog.vao = 0;
    assert(check_glerror());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_gluniform_location(const glprogram& prog, const std::string& name) {
    assert(check_glerror());
    return glGetUniformLocation(prog.pid, name.c_str());
    assert(check_glerror());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_glattrib_location(const glprogram& prog, const std::string& name) {
    assert(check_glerror());
    return glGetAttribLocation(prog.pid, name.c_str());
    assert(check_glerror());
}

// Set uniform values.
void set_gluniform(const glprogram& prog, int pos, bool val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1i(pos, (int)val);
}
void set_gluniform(const glprogram& prog, int pos, int val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1i(pos, val);
}
void set_gluniform(const glprogram& prog, int pos, float val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1f(pos, val);
}
void set_gluniform(const glprogram& prog, int pos, const vec2f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform2f(pos, val.x, val.y);
}
void set_gluniform(const glprogram& prog, int pos, const vec3f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform3f(pos, val.x, val.y, val.z);
}
void set_gluniform(const glprogram& prog, int pos, const vec4f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform4f(pos, val.x, val.y, val.z, val.w);
}
void set_gluniform(const glprogram& prog, int pos, const mat4f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniformMatrix4fv(pos, 1, false, &val.x.x);
}
void set_gluniform(const glprogram& prog, int pos, const frame3f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniformMatrix4x3fv(pos, 1, false, &val.x.x);
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture_info& tinfo, uint tunit) {
    static const auto wrap_mode_map = std::unordered_map<gltexture_wrap, uint>{
        {gltexture_wrap::repeat, GL_REPEAT},
        {gltexture_wrap::clamp, GL_CLAMP_TO_EDGE},
        {gltexture_wrap::mirror, GL_MIRRORED_REPEAT}};
    static const auto filter_mode_map =
        std::unordered_map<gltexture_filter, uint>{
            {gltexture_filter::nearest, GL_NEAREST},
            {gltexture_filter::linear, GL_LINEAR},
            {gltexture_filter::nearest_mipmap_nearest,
                GL_NEAREST_MIPMAP_NEAREST},
            {gltexture_filter::linear_mipmap_nearest, GL_LINEAR_MIPMAP_NEAREST},
            {gltexture_filter::nearest_mipmap_linear, GL_NEAREST_MIPMAP_LINEAR},
            {gltexture_filter::linear_mipmap_linear, GL_LINEAR_MIPMAP_LINEAR}};

    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_gltexture_valid(tinfo.txt)) {
        bind_gltexture(tinfo.txt, tunit);
        if (tinfo.wrap_s != gltexture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                wrap_mode_map.at(tinfo.wrap_s));
        if (tinfo.wrap_t != gltexture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                wrap_mode_map.at(tinfo.wrap_t));
        if (tinfo.filter_min != gltexture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                filter_mode_map.at(tinfo.filter_min));
        if (tinfo.filter_mag != gltexture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                filter_mode_map.at(tinfo.filter_mag));
        glUniform1i(pos, tunit);
    } else {
        unbind_gltexture(tinfo.txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(check_glerror());
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture& txt, uint tunit) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_gltexture_valid(txt)) {
        bind_gltexture(txt, tunit);
        glUniform1i(pos, tunit);
    } else {
        unbind_gltexture(txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(check_glerror());
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(const glprogram& prog, const std::string& name,
    const gltexture& txt, uint tunit) {
    set_gluniform_texture(prog, get_gluniform_location(prog, name), txt, tunit);
}

// Sets a vartex attribute for program pid and variable var.
void set_glattribute(
    const glprogram& prog, int pos, const glvertex_buffer& buf, float def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib1f(pos, def);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec2f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib2f(pos, def.x, def.y);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec3f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib3f(pos, def.x, def.y, def.z);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec4f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib4f(pos, def.x, def.y, def.z, def.w);
    }
    assert(check_glerror());
}

// Binds a program
void bind_glprogram(const glprogram& prog) {
    assert(check_glerror());
    if (!prog.pid) return;
    glBindVertexArray(prog.vao);
    glUseProgram(prog.pid);
    assert(check_glerror());
}

// Unbind a program
void unbind_glprogram(const glprogram& prog) {
    assert(check_glerror());
    glUseProgram(0);
    glBindVertexArray(0);
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

static std::string glimage_vert =
    R"(
    #version 330

    layout(location = 0) in vec2 vert_texcoord;

    uniform vec2 offset;
    uniform vec2 win_size;
    uniform float zoom;
    uniform sampler2D img;

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

static std::string glimage_frag =
    R"(
    #version 330

    in vec2 texcoord;

    uniform sampler2D img;
    uniform float exposure;
    uniform float gamma;

    out vec4 color;

    void main() {
        vec4 c = texture(img,texcoord);
        c.xyz = pow(c.xyz * pow(2,exposure), vec3(1/gamma));
        color = c;
    }
    )";

// Initialize the program. Call with true only after the GL is initialized.
glimage_program make_glimage_program() {
    auto prog = glimage_program();
    prog.prog = make_glprogram(glimage_vert, glimage_frag);

    update_glbuffer(
        prog.vbo, false, std::vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}});
    update_glbuffer(prog.ebo, true, std::vector<vec3i>{{0, 1, 2}, {0, 2, 3}});
    return prog;
}

// Draws the stdimage program.
void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma) {
    assert(is_gltexture_valid(txt));

    bind_glprogram(prog.prog);

    enable_glblending(true);
    set_glblend_over();

    bind_gltexture(txt, 0);
    set_gluniform(prog.prog, "zoom", zoom);
    set_gluniform(
        prog.prog, "win_size", vec2f{(float)win_size.x, (float)win_size.y});
    set_gluniform(prog.prog, "offset", offset);
    set_gluniform(prog.prog, "exposure", exposure);
    set_gluniform(prog.prog, "gamma", gamma);
    set_gluniform_texture(prog.prog, "img", txt, 0);

    set_glattribute(prog.prog, "vert_texcoord", prog.vbo, vec2f{0, 0});
    draw_glelems(prog.ebo);

    unbind_glprogram(prog.prog);

    enable_glblending(false);

    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL STANDARD SURFACE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize a standard shader. Call with true only after the gl has
// been initialized
glsurface_program make_glsurface_program() {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    std::string _vert =
        R"(
        #version 330

        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;           // shape transform
        uniform float shape_normal_offset;           // shape normal offset

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

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

            // normal offset
            if(shape_normal_offset != 0) {
                pos += shape_normal_offset * norm;
            }

            // world projection
            pos = (shape_xform * vec4(pos,1)).xyz;
            norm = (shape_xform * vec4(norm,0)).xyz;
            tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
        }
        )";

    std::string _frag_header =
        R"(
        #version 330

        float pi = 3.14159265;

        )";

    std::string _frag_lighting =
        R"(
        uniform bool eyelight;        // eyelight shading
        uniform vec3 lamb;             // ambient light
        uniform int lnum;              // number of lights
        uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
        uniform vec3 lpos[16];         // light positions
        uniform vec3 lke[16];          // light intensities

        void eval_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(ltype[lid] == 0) {
                // compute point light color at pos
                cl = lke[lid] / pow(length(lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lpos[lid]-pos);
            }
            else if(ltype[lid] == 1) {
                // compute light color
                cl = lke[lid];
                // compute light direction
                wi = normalize(lpos[lid]);
            }
        }

        )";

    std::string _frag_brdf =
        R"(
        vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
            vec3 n, vec3 wi, vec3 wo) {
            if(etype == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(rs*rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(etype == 1) {
                return ((1+dot(wo,wi))/2) * kd/pi;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pi;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pi)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pi;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pi)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        )";

    std::string _frag_material =
        R"(
        uniform int elem_type;
        uniform bool elem_faceted;
        uniform vec4 highlight;   // highlighted color

        uniform int mat_type;          // material type
        uniform vec3 mat_ke;           // material ke
        uniform vec3 mat_kd;           // material kd
        uniform vec3 mat_ks;           // material ks
        uniform float mat_rs;          // material rs
        uniform float mat_op;          // material op

        uniform bool mat_txt_ke_on;    // material ke texture on
        uniform sampler2D mat_txt_ke;  // material ke texture
        uniform bool mat_txt_kd_on;    // material kd texture on
        uniform sampler2D mat_txt_kd;  // material kd texture
        uniform bool mat_txt_ks_on;    // material ks texture on
        uniform sampler2D mat_txt_ks;  // material ks texture
        uniform bool mat_txt_rs_on;    // material rs texture on
        uniform sampler2D mat_txt_rs;  // material rs texture

        uniform bool mat_txt_norm_on;    // material norm texture on
        uniform sampler2D mat_txt_norm;  // material norm texture
        uniform float mat_txt_norm_scale;  // material norm scale

        uniform bool mat_txt_occ_on;    // material occ texture on
        uniform sampler2D mat_txt_occ;  // material occ texture
        uniform float mat_txt_occ_scale;  // material occ scale

        uniform bool mat_use_phong;       // material use phong
        uniform bool mat_double_sided;    // material double sided
        uniform bool mat_alpha_cutout;    // material alpha cutout

        bool eval_material(vec2 texcoord, vec4 color, out vec3 ke, 
                           out vec3 kd, out vec3 ks, out float rs, out float op, out bool cutout) {
            if(mat_type == 0) {
                ke = mat_ke;
                kd = vec3(0,0,0);
                ks = vec3(0,0,0);
                op = 1;
                return false;
            }

            ke = color.xyz * mat_ke;
            kd = color.xyz * mat_kd;
            ks = color.xyz * mat_ks;
            rs = mat_rs;
            op = color.w * mat_op;

            vec3 ke_txt = (mat_txt_ke_on) ? texture(mat_txt_ke,texcoord).xyz : vec3(1);
            vec4 kd_txt = (mat_txt_kd_on) ? texture(mat_txt_kd,texcoord) : vec4(1);
            vec4 ks_txt = (mat_txt_ks_on) ? texture(mat_txt_ks,texcoord) : vec4(1);
            float rs_txt = (mat_txt_rs_on) ? texture(mat_txt_rs,texcoord).x : 1;

            // scale common values
            ke *= ke_txt;

            // get material color from textures and adjust values
            if(mat_type == 0) {
            } else if(mat_type == 1) {
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt;
                rs = rs*rs;
            } else if(mat_type == 2) {
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(mat_type == 3) {
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= ks_txt.w;
                rs = (1 - rs) * (1 - rs);
                op *= kd_txt.w;
            }

            cutout = mat_alpha_cutout && op == 0;
            return true;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!mat_txt_norm_on) return norm;
            vec3 tangu = normalize(tangsp.xyz);
            vec3 tangv = normalize(cross(tangu, norm));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 txt = 2 * pow(texture(mat_txt_norm,texcoord).xyz, vec3(1/2.2)) - 1;
            return normalize( tangu * txt.x + tangv * txt.y + norm * txt.z );
        }

        )";

    std::string _frag_main =
        R"(
        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        uniform float exposure; 
        uniform float gamma;

        out vec4 frag_color;      

        vec3 triangle_normal(vec3 pos) {
            vec3 fdx = dFdx(pos); 
            vec3 fdy = dFdy(pos); 
            return normalize(cross(fdx, fdy));
        }

        // main
        void main() {
            // view vector
            vec3 wo = normalize( (cam_xform*vec4(0,0,0,1)).xyz - pos );

            // prepare normals
            vec3 n;
            if(elem_faceted) {
                n = triangle_normal(pos);
            } else {
                n = normalize(norm);
            }

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(mat_double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op; bool brdf_cutout;
            bool has_brdf = eval_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, brdf_cutout);

            // exit if needed
            if(brdf_cutout) discard;

            // check const color
            if(elem_type == 0) {
                frag_color = vec4(brdf_ke,brdf_op);
                return;
            }

            // emission
            vec3 c = brdf_ke;

            // check early exit
            if(brdf_kd != vec3(0,0,0) || brdf_ks != vec3(0,0,0)) {
                // eyelight shading
                if(eyelight) {
                    vec3 wi = wo;
                    c += pi * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lamb * brdf_kd;
                    // foreach light
                    for(int lid = 0; lid < lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        eval_light(lid, pos, cl, wi);
                        c += cl * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                    }
                }
            }

            // final color correction
            c = pow(c * pow(2,exposure), vec3(1/gamma));

            // highlighting
            if(highlight.w > 0) {
                if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
                    c = highlight.xyz * highlight.w + c * (1-highlight.w);
            }

            // output final color by setting gl_FragColor
            frag_color = vec4(c,brdf_op);
        }
        )";
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

    assert(check_glerror());
    auto prog = glsurface_program();
    prog.prog =
        make_glprogram(_vert, _frag_header + _frag_lighting + _frag_brdf +
                                  _frag_material + _frag_main);
    assert(check_glerror());
    prog.eyelight_id = get_gluniform_location(prog.prog, "eyelight");
    prog.exposure_id = get_gluniform_location(prog.prog, "exposure");
    prog.gamma_id = get_gluniform_location(prog.prog, "gamma");
    prog.cam_xform_id = get_gluniform_location(prog.prog, "cam_xform");
    prog.cam_xform_inv_id = get_gluniform_location(prog.prog, "cam_xform_inv");
    prog.cam_proj_id = get_gluniform_location(prog.prog, "cam_proj");
    prog.lamb_id = get_gluniform_location(prog.prog, "lamb");
    prog.lnum_id = get_gluniform_location(prog.prog, "lnum");
    for (auto i = 0; i < 16; i++) {
        auto is = std::to_string(i);
        prog.lpos_id[i] = get_gluniform_location(prog.prog, "lpos[" + is + "]");
        prog.lke_id[i] = get_gluniform_location(prog.prog, "lke[" + is + "]");
        prog.ltype_id[i] =
            get_gluniform_location(prog.prog, "ltype[" + is + "]");
    }
    prog.shp_xform_id = get_gluniform_location(prog.prog, "shape_xform");
    prog.shp_normal_offset_id =
        get_gluniform_location(prog.prog, "shape_normal_offset");
    prog.highlight_id = get_gluniform_location(prog.prog, "highlight");
    prog.mtype_id = get_gluniform_location(prog.prog, "mat_type");
    prog.ke_id = get_gluniform_location(prog.prog, "mat_ke");
    prog.kd_id = get_gluniform_location(prog.prog, "mat_kd");
    prog.ks_id = get_gluniform_location(prog.prog, "mat_ks");
    prog.rs_id = get_gluniform_location(prog.prog, "mat_rs");
    prog.op_id = get_gluniform_location(prog.prog, "mat_op");
    prog.ke_txt_id = get_gluniform_location(prog.prog, "mat_txt_ke");
    prog.ke_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_ke_on");
    prog.kd_txt_id = get_gluniform_location(prog.prog, "mat_txt_kd");
    prog.kd_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_kd_on");
    prog.ks_txt_id = get_gluniform_location(prog.prog, "mat_txt_ks");
    prog.ks_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_ks_on");
    prog.rs_txt_id = get_gluniform_location(prog.prog, "mat_txt_rs");
    prog.rs_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_rs_on");
    prog.norm_txt_id = get_gluniform_location(prog.prog, "mat_txt_norm");
    prog.norm_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_norm_on");
    prog.occ_txt_id = get_gluniform_location(prog.prog, "mat_txt_occ");
    prog.occ_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_occ_on");
    prog.norm_scale_id = get_gluniform_location(prog.prog, "mat_norm_scale");
    prog.occ_scale_id = get_gluniform_location(prog.prog, "mat_occ_scale");
    prog.use_phong_id = get_gluniform_location(prog.prog, "mat_use_phong");
    prog.double_sided_id =
        get_gluniform_location(prog.prog, "mat_double_sided");
    prog.alpha_cutout_id =
        get_gluniform_location(prog.prog, "mat_alpha_cutout");
    prog.etype_id = get_gluniform_location(prog.prog, "elem_type");
    prog.efaceted_id = get_gluniform_location(prog.prog, "elem_faceted");
    assert(check_glerror());
    prog.pos_id = get_glattrib_location(prog.prog, "vert_pos");
    prog.norm_id = get_glattrib_location(prog.prog, "vert_norm");
    prog.texcoord_id = get_glattrib_location(prog.prog, "vert_texcoord");
    prog.color_id = get_glattrib_location(prog.prog, "vert_color");
    prog.tangsp_id = get_glattrib_location(prog.prog, "vert_tangsp");
    assert(check_glerror());
    return prog;
}

// Starts a frame by setting exposure/gamma values, camera transforms
// and projection. Sets also whether to use full shading or a quick
// eyelight preview.
void begin_glsurface_frame(const glsurface_program& prog,
    const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj, bool shade_eyelight, float exposure,
    float gamma) {
    assert(check_glerror());
    bind_glprogram(prog.prog);
    set_gluniform(prog.prog, prog.cam_xform_id, camera_xform);
    set_gluniform(prog.prog, prog.cam_xform_inv_id, camera_xform_inv);
    set_gluniform(prog.prog, prog.cam_proj_id, camera_proj);
    set_gluniform(prog.prog, prog.eyelight_id, shade_eyelight);
    set_gluniform(prog.prog, prog.exposure_id, exposure);
    set_gluniform(prog.prog, prog.gamma_id, gamma);
    assert(check_glerror());
    set_glsurface_elems(prog, glelem_type::triangle, false);
}

// Ends a frame.
void end_glsurface_frame(const glsurface_program& prog) {
    assert(check_glerror());
    unbind_glprogram(prog.prog);
    //    glBindVertexArray(0);
    //    glUseProgram(0);
    assert(check_glerror());
}

// Set num lights with position pos, color ke, type ltype. Also set the
// ambient illumination amb.
void set_glsurface_lights(
    const glsurface_program& prog, const vec3f& amb, const gllights& lights) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.lamb_id, amb);
    set_gluniform(prog.prog, prog.lnum_id, (int)lights.pos.size());
    for (auto i = 0; i < lights.pos.size(); i++) {
        set_gluniform(prog.prog, prog.lpos_id[i], lights.pos[i]);
        set_gluniform(prog.prog, prog.lke_id[i], lights.ke[i]);
        set_gluniform(prog.prog, prog.ltype_id[i], (int)lights.type[i]);
    }
    assert(check_glerror());
}

// Begins drawing a shape with transform xform.
void begin_glsurface_shape(
    const glsurface_program& prog, const mat4f& xform, float normal_offset) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.shp_xform_id, xform);
    set_gluniform(prog.prog, prog.shp_normal_offset_id, normal_offset);
    assert(check_glerror());
}

// End shade drawing.
void end_glsurface_shape(const glsurface_program& prog) {
    assert(check_glerror());
    for (int i = 0; i < 16; i++) unbind_glbuffer(i);
    assert(check_glerror());
}

// Sets normal offset.
void set_glsurface_normaloffset(
    const glsurface_program& prog, float normal_offset) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.shp_normal_offset_id, normal_offset);
    assert(check_glerror());
}

// Set the object as highlighted.
void set_glsurface_highlight(
    const glsurface_program& prog, const vec4f& highlight) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.highlight_id, highlight);
    assert(check_glerror());
}

// Set material values with emission ke, diffuse kd, specular ks and
// specular roughness rs, opacity op. Indicates textures ids with the
// correspoinding XXX_txt variables. Sets also normal and occlusion
// maps. Works for points/lines/triangles (diffuse for points,
// Kajiya-Kay for lines, GGX/Phong for triangles).
// Material type matches the scene material type.
void set_glsurface_material(const glsurface_program& prog, int mtype,
    const vec3f& ke, const vec3f& kd, const vec3f& ks, float rs, float op,
    const gltexture_info& ke_txt, const gltexture_info& kd_txt,
    const gltexture_info& ks_txt, const gltexture_info& rs_txt,
    const gltexture_info& norm_txt, const gltexture_info& occ_txt,
    bool use_phong, bool double_sided, bool alpha_cutout) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.mtype_id, mtype);
    set_gluniform(prog.prog, prog.ke_id, ke);
    set_gluniform(prog.prog, prog.kd_id, kd);
    set_gluniform(prog.prog, prog.ks_id, ks);
    set_gluniform(prog.prog, prog.rs_id, rs);
    set_gluniform(prog.prog, prog.op_id, op);
    set_gluniform_texture(
        prog.prog, prog.ke_txt_id, prog.ke_txt_on_id, ke_txt, 0);
    set_gluniform_texture(
        prog.prog, prog.kd_txt_id, prog.kd_txt_on_id, kd_txt, 1);
    set_gluniform_texture(
        prog.prog, prog.ks_txt_id, prog.ks_txt_on_id, ks_txt, 2);
    set_gluniform_texture(
        prog.prog, prog.rs_txt_id, prog.rs_txt_on_id, rs_txt, 3);
    set_gluniform_texture(
        prog.prog, prog.norm_txt_id, prog.norm_txt_on_id, norm_txt, 4);
    //    set_gluniform(prog.prog, prog.norm_scale_id, norm_txt.scale);
    //    set_gluniform(prog.prog, prog.occ_scale_id, occ_txt.scale);
    //    set_gluniform(prog.prog, prog.use_phong_id, use_phong);
    set_gluniform(prog.prog, prog.double_sided_id, double_sided);
    set_gluniform(prog.prog, prog.alpha_cutout_id, alpha_cutout);
    assert(check_glerror());
}

// Set constant material values with emission ke.
void set_glsurface_constmaterial(
    const glsurface_program& prog, const vec3f& ke, float op) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.mtype_id, 0);
    set_gluniform(prog.prog, prog.ke_id, ke);
    set_gluniform(prog.prog, prog.op_id, op);
    assert(check_glerror());
}

// Set vertex data with buffers for position pos, normals norm, texture
// coordinates texcoord, per-vertex color color and tangent space
// tangsp.
void set_glsurface_vert(const glsurface_program& prog,
    const glvertex_buffer& pos, const glvertex_buffer& norm,
    const glvertex_buffer& texcoord, const glvertex_buffer& color,
    const glvertex_buffer& tangsp) {
    assert(check_glerror());
    set_glattribute(prog.prog, prog.pos_id, pos, zero3f);
    set_glattribute(prog.prog, prog.norm_id, norm, zero3f);
    set_glattribute(prog.prog, prog.texcoord_id, texcoord, zero2f);
    set_glattribute(prog.prog, prog.color_id, color, vec4f{1, 1, 1, 1});
    set_glattribute(prog.prog, prog.tangsp_id, tangsp, zero4f);
    assert(check_glerror());
}

// Set element properties.
void set_glsurface_elems(
    const glsurface_program& prog, glelem_type etype, bool faceted) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.etype_id, (int)etype);
    set_gluniform(prog.prog, prog.efaceted_id, (int)faceted);
    assert(check_glerror());
}

glwindow::~glwindow() {
    if (widget_enabled) {
        ImGui_ImplGlfwGL3_Shutdown();
        widget_enabled = false;
    }
    if (gwin) {
        glfwDestroyWindow(gwin);
        glfwTerminate();
        gwin = nullptr;
    }
}

// Support
void _glfw_error_cb(int error, const char* description) {
    log_error("GLFW error: {}\n", description);
}

// Support
void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) { ImGui_ImplGlfw_CharCallback(win->gwin, key); }
    if (win->text_cb) win->text_cb(key);
}

// Support
void _glfw_key_cb(
    GLFWwindow* gwin, int key, int scancode, int action, int mods) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_KeyCallback(win->gwin, key, scancode, action, mods);
    }
}

// Support
void _glfw_mouse_cb(GLFWwindow* gwin, int button, int action, int mods) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_MouseButtonCallback(win->gwin, button, action, mods);
    }
    if (win->mouse_cb) win->mouse_cb(button, action == GLFW_PRESS, mods);
}

// Support
void _glfw_scroll_cb(GLFWwindow* gwin, double xoffset, double yoffset) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_ScrollCallback(win->gwin, xoffset, yoffset);
    }
}

// Support
void _glfw_refresh_cb(GLFWwindow* gwin) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->refresh_cb) win->refresh_cb();
}

// Initialize glwindow
glwindow* make_glwindow(
    int width, int height, const std::string& title, bool opengl4) {
    auto win = new glwindow();

    // glwindow
    glfwSetErrorCallback(_glfw_error_cb);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    // profile creation
    if (opengl4) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    } else {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    }
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    win->gwin = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->gwin);
    glfwSetWindowUserPointer(win->gwin, win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win->gwin, _glfw_text_cb);
    glfwSetKeyCallback(win->gwin, _glfw_key_cb);
    glfwSetMouseButtonCallback(win->gwin, _glfw_mouse_cb);
    glfwSetScrollCallback(win->gwin, _glfw_scroll_cb);

    glfwSetWindowRefreshCallback(win->gwin, _glfw_refresh_cb);

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif
    return win;
}

// Set glwindow callbacks
void set_glwindow_callbacks(glwindow* win, text_glcallback text_cb,
    mouse_glcallback mouse_cb, refresh_glcallback refresh_cb) {
    win->text_cb = text_cb;
    win->mouse_cb = mouse_cb;
    win->refresh_cb = refresh_cb;
    if (win->text_cb) glfwSetCharCallback(win->gwin, _glfw_text_cb);
}

// Set glwindow title
void set_glwindow_title(glwindow* win, const std::string& title) {
    glfwSetWindowTitle(win->gwin, title.c_str());
}

// Wait events
void wait_glwindow_events(glwindow* win) { glfwWaitEvents(); }

// Poll events
void poll_glwindow_events(glwindow* win) { glfwPollEvents(); }

#ifdef __APPLE__

// Wait events
void wait_glwindow_events_timeout(glwindow* win, double timeout_sec) {
    glfwWaitEventsTimeout(timeout_sec);
}

// Wait events
void post_glwindow_event(glwindow* win) { glfwPostEmptyEvent(); }

#endif

// Swap buffers
void swap_glwindow_buffers(glwindow* win) { glfwSwapBuffers(win->gwin); }

// Should close
bool should_glwindow_close(glwindow* win) {
    return glfwWindowShouldClose(win->gwin);
}

// Mouse button
int get_glmouse_button(glwindow* win) {
    auto mouse1 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
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

// Mouse position
vec2i get_glmouse_pos(glwindow* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(int)x, (int)y};
}

// Mouse position
vec2f get_glmouse_posf(glwindow* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(float)x, (float)y};
}

// Window size
vec2i get_glwindow_size(glwindow* win) {
    auto ret = vec2i{0, 0};
    glfwGetWindowSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Check if a key is pressed (not all keys are supported)
bool get_glkey(glwindow* win, int key) {
    key = std::toupper(key);
    return glfwGetKey(win->gwin, key) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glalt_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glctrl_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glshift_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

// Framebuffer size
vec2i get_glframebuffer_size(glwindow* win) {
    auto ret = vec2i{0, 0};
    glfwGetFramebufferSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Read pixels
void take_glscreenshot4b(glwindow* win, int& width, int& height,
    std::vector<ygl::vec4b>& img, bool flipy, bool back) {
    auto wh = get_glframebuffer_size(win);
    width = wh.x;
    height = wh.y;
    img = std::vector<ygl::vec4b>(wh.x * wh.y);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    if (flipy) {
        for (int j = 0; j < height / 2; j++) {
            for (auto i = 0; i < width; i++) {
                std::swap(
                    img[i + width * j], img[i + width * (height - 1 - j)]);
            }
        }
    }
}

// Initialize widgets
void init_imgui(glwindow* win, bool light_style, bool alt_font) {
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(win->gwin, false);
    ImGuiIO& io = ImGui::GetIO();
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard
    // Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable
    // Gamepad Controls
    ImGui::GetStyle().WindowRounding = 0;
    io.IniFilename = nullptr;
    auto size = get_glwindow_size(win);
    ImGui::SetNextWindowPos({(float)size.x - 320, 0});
    ImGui::SetNextWindowSize({(float)320, (float)size.y});
    if (light_style) ImGui::StyleColorsLight();
    if (alt_font) {
        io.Fonts->AddFontFromMemoryCompressedTTF(
            imgui_extrafont_compressed_data(),
            imgui_extrafont_compressed_size(), 16);
    } else {
        io.Fonts->AddFontDefault();
    }
    win->widget_enabled = true;
}

// Begin draw widget
bool begin_imgui_frame(glwindow* win, const std::string& title) {
    static bool first_time = true;
    ImGui_ImplGlfwGL3_NewFrame();
    if (first_time) {
        auto size = get_glwindow_size(win);
        ImGui::SetNextWindowPos({(float)size.x - 320, 0});
        ImGui::SetNextWindowSize({(float)320, (float)size.y});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    ImGui::Begin(title.c_str(), nullptr);
    // ImGui::ShowTestWindow();
    // ImGui::ShowStyleEditor();
    return true;
}

// End draw widget
void end_imgui_frame(glwindow* win) {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
}

// Whether widget are active
bool get_imgui_active(glwindow* win) {
    if (!win->widget_enabled) return false;
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

// Horizontal separator
void draw_imgui_separator(glwindow* win) { ImGui::Separator(); }

// Indent widget
void begin_imgui_indent(glwindow* win) { ImGui::Indent(); }

// Indent widget
void end_imgui_indent(glwindow* win) { ImGui::Unindent(); }

// Continue line with next widget
void continue_imgui_line(glwindow* win) { ImGui::SameLine(); }

// Label widget
void draw_imgui_label(
    glwindow* win, const std::string& lbl, const std::string& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg.c_str());
}

// Value widget
bool draw_imgui_text(glwindow* win, const std::string& lbl, std::string& str) {
    char buf[4096];
    if (str.length() >= 4096) throw std::runtime_error("bad memory");
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputText(lbl.c_str(), buf, 4096);
    str = buf;
    return ret;
}

// Value widget
bool draw_imgui_multiline_text(
    glwindow* win, const std::string& lbl, std::string& str) {
    char sbuf[8192];
    std::vector<char> dbuf;
    char* buf = nullptr;
    int buf_size = 0;
    if (str.size() > sizeof(sbuf) / 2) {
        dbuf.resize(str.size() * 2);
        buf = dbuf.data();
        buf_size = dbuf.size();
    } else {
        buf = sbuf;
        buf_size = sizeof(sbuf);
    }
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputTextMultiline(lbl.c_str(), buf, buf_size);
    str = buf;
    return ret;
}

static float draw_drag_scale = 1 / 100.0f;

// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, int& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt(lbl.c_str(), &val, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec2i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt2(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec3i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt3(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec4i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt4(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, float& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat(lbl.c_str(), &val, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec2f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat2(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec3f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat3(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec4f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat4(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, mat4f& val, float min, float max) {
    auto modx = draw_imgui_dragbox(win, lbl + ".x", val.x, min, max);
    auto mody = draw_imgui_dragbox(win, lbl + ".y", val.y, min, max);
    auto modz = draw_imgui_dragbox(win, lbl + ".z", val.z, min, max);
    auto modw = draw_imgui_dragbox(win, lbl + ".w", val.w, min, max);
    return modx || mody || modz || modw;
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, frame3f& val, float min, float max) {
    auto modx = draw_imgui_dragbox(win, lbl + ".x", val.x, -1, 1);
    auto mody = draw_imgui_dragbox(win, lbl + ".y", val.y, -1, 1);
    auto modz = draw_imgui_dragbox(win, lbl + ".z", val.z, -1, 1);
    auto modo = draw_imgui_dragbox(win, lbl + ".o", val.o, min, max);
    return modx || mody || modz || modo;
}

// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec3f& val) {
    auto mod = ImGui::ColorEdit3(
        lbl.c_str(), (float*)&val.x, ImGuiColorEditFlags_Float);
    // fix for bug in ImGui
    if (mod) {
        if (val.x < 0.0001f) val.x = 0;
        if (val.y < 0.0001f) val.y = 0;
        if (val.z < 0.0001f) val.z = 0;
    }
    return mod;
}
// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4f& val) {
    auto mod = ImGui::ColorEdit4(
        lbl.c_str(), (float*)&val.x, ImGuiColorEditFlags_Float);
    // fix for bug in ImGui
    if (mod) {
        if (val.x < 0.0001f) val.x = 0;
        if (val.y < 0.0001f) val.y = 0;
        if (val.z < 0.0001f) val.z = 0;
        if (val.w < 0.0001f) val.w = 0;
    }
    return mod;
}
// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4b& val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)&val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)&val = valb;
        return true;
    }
    return false;
}
bool draw_hdr_color_widget(
    glwindow* win, const std::string& lbl, vec3f& val, float max) {
    auto vall = ygl::max(val);
    auto valc = val;
    if (vall > 1) {
        valc /= vall;
    } else {
        vall = 1;
    }
    auto mod1 = draw_imgui_dragbox(win, lbl + "(m)", vall, 0, max);
    auto mod2 = draw_imgui_colorbox(win, lbl, valc);
    if (mod1 || mod2) {
        val = valc * vall;
        return true;
    } else {
        return false;
    }
}

// Bool widget
bool draw_imgui_checkbox(glwindow* win, const std::string& lbl, bool& val) {
    return ImGui::Checkbox(lbl.c_str(), &val);
}

// Combo widget
bool begin_imgui_combobox(
    glwindow* win, const std::string& lbl, const std::string& label) {
    return ImGui::BeginCombo(lbl.c_str(), label.c_str());
}

// Combo widget
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, bool selected) {
    ImGui::PushID((void*)(intptr_t)idx);
    auto clicked = ImGui::Selectable(label.c_str(), selected);
    if (selected) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
    return clicked;
}

// Combo widget
void end_imgui_combobox(glwindow* win) { ImGui::EndCombo(); }

// Button widget
bool draw_imgui_button(glwindow* win, const std::string& lbl) {
    return ImGui::Button(lbl.c_str());
}

// Collapsible header
bool draw_imgui_header(glwindow* win, const std::string& lbl) {
    return ImGui::CollapsingHeader(lbl.c_str());
}

// Start tree node
bool begin_imgui_tree(glwindow* win, const std::string& lbl) {
    return ImGui::TreeNode(lbl.c_str());
}

// Collapsible header
void end_imgui_tree(glwindow* win) { ImGui::TreePop(); }

// Start selectable tree node
bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
    return open;
}

// End selectable tree node
void end_imgui_tree(glwindow* win, void* content) { ImGui::TreePop(); }

// Selectable tree leaf node
void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
}

// Selectable tree leaf node
void draw_imgui_tree_leaf(glwindow* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    draw_imgui_tree_leaf(win, lbl, selection, content);
    ImGui::PopStyleColor();
}

// Image widget
void draw_imgui_imagebox(
    glwindow* win, int tid, const vec2i& size, const vec2i& imsize) {
    auto w = ImGui::GetContentRegionAvailWidth();
    auto s = vec2f{(float)size.x, (float)size.y};
    auto a = (float)imsize.x / (float)imsize.y;
    if (!s.x && !s.y) {
        s.x = w;
        s.y = w / a;
    } else if (s.x && !s.y) {
        s.y = s.x / a;
    } else if (!s.x && s.y) {
        s.x = s.y * a;
    } else {
        auto as = s.x / s.y;
        if (as / a > 1) {
            s.x = s.y * a;
        } else {
            s.y = s.x / a;
        }
    }
    if (s.x > w) {
        s.x = w;
        s.y = w / a;
    }
    ImGui::Image((void*)(size_t)tid, {s.x, s.y});
}

// Image widget
void draw_imgui_imagebox(
    glwindow* win, const gltexture& txt, const vec2i& size) {
    draw_imgui_imagebox(
        win, get_gltexture_id(txt), size, {txt.width, txt.height});
}

// Scroll region
void begin_imgui_scrollarea(
    glwindow* win, const std::string& lbl, int height, bool border) {
    ImGui::BeginChild(lbl.c_str(), ImVec2(0, height), border);
}
// Scroll region
void end_imgui_scrollarea(glwindow* win) { ImGui::EndChild(); }
// Scroll region
void move_imgui_scrollarea(glwindow* win) { ImGui::SetScrollHere(); }

// Group ids
void push_imgui_groupid(glwindow* win, int gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, void* gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, const void* gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, const char* gid) { ImGui::PushID(gid); }
// Group ids
void pop_imgui_groupid(glwindow* win) { ImGui::PopID(); }

// Widget style
void push_imgui_style(glwindow* win, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
}

// Widget style
void pop_imgui_style(glwindow* win) { ImGui::PopStyleColor(); }

}  // namespace ygl
