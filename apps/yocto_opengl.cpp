//
// Utilities to use OpenGL 3, GLFW and ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_opengl.h"

#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"

#include <algorithm>
#include <cstdarg>
#include <deque>
#include <mutex>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"

#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"

#define CUTE_FILES_IMPLEMENTATION
#include "ext/cute_files.h"

// -----------------------------------------------------------------------------
// LOW-LEVEL OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

void check_glerror() {
  if (glGetError() != GL_NO_ERROR) printf("gl error\n");
}

void clear_glframebuffer(const vec4f& color, bool clear_depth) {
  glClearColor(color.x, color.y, color.z, color.w);
  if (clear_depth) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
  } else {
    glClear(GL_COLOR_BUFFER_BIT);
  }
}

void set_glviewport(const vec4i& viewport) {
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);
}

void set_glwireframe(bool enabled) {
  if (enabled)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void set_glblending(bool enabled) {
  if (enabled) {
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
  } else {
    glDisable(GL_BLEND);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// init image program
void init_glimage_program(opengl_image& glimage) {
  if (glimage.program) return;
  auto vert =
      R"(
      #version 330
      in vec2 texcoord;
      out vec2 frag_texcoord;
      uniform vec2 window_size, image_size;
      uniform vec2 image_center;
      uniform float image_scale;
      void main() {
          vec2 pos = (texcoord - vec2(0.5,0.5)) * image_size * image_scale + image_center;
          gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
          frag_texcoord = texcoord;
      }
      )";
#if 0
  auto vert = R"(
          #version 330
          in vec2 texcoord;
          out vec2 frag_texcoord;
          uniform vec2 window_size, image_size, border_size;
          uniform vec2 image_center;
          uniform float image_scale;
          void main() {
              vec2 pos = (texcoord - vec2(0.5,0.5)) * (image_size + border_size*2) * image_scale + image_center;
              gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
              frag_texcoord = texcoord;
          }
      )";
#endif
  auto frag =
      R"(
      #version 330
      in vec2 frag_texcoord;
      out vec4 frag_color;
      uniform sampler2D txt;
      void main() {
          frag_color = texture(txt, frag_texcoord);
      }
      )";
#if 0
    auto frag = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform vec2 image_size, border_size;
            uniform float image_scale;
            void main() {
                ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
                ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
                ivec2 tile = tilecoord / 16;
                if(imcoord.x <= 0 || imcoord.y <= 0 || 
                    imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
                else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
                else frag_color = vec4(0.3,0.3,0.3,1);
            }
        )";
#endif

  init_glprogram(glimage.program, vert, frag);
  init_glarraybuffer(
      glimage.texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
  init_glelementbuffer(
      glimage.element, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
}

// update image data
void update_glimage(
    opengl_image& glimage, const image<vec4f>& img, bool linear, bool mipmap) {
  init_glimage_program(glimage);
  if (!glimage.texture) {
    init_gltexture(glimage.texture, img, false, linear, mipmap);
  } else if (glimage.texture.size != img.size() ||
             glimage.texture.linear != linear ||
             glimage.texture.mipmap != mipmap) {
    delete_gltexture(glimage.texture);
    init_gltexture(glimage.texture, img, false, linear, mipmap);
  } else {
    update_gltexture(glimage.texture, img, mipmap);
  }
}
void update_glimage(
    opengl_image& glimage, const image<vec4b>& img, bool linear, bool mipmap) {
  init_glimage_program(glimage);
  if (!glimage.texture) {
    init_gltexture(glimage.texture, img, false, linear, mipmap);
  } else if (glimage.texture.size != img.size() ||
             glimage.texture.linear != linear ||
             glimage.texture.mipmap != mipmap) {
    delete_gltexture(glimage.texture);
    init_gltexture(glimage.texture, img, false, linear, mipmap);
  } else {
    update_gltexture(glimage.texture, img, mipmap);
  }
}

void update_glimage_region(opengl_image& glimage, const image<vec4f>& img,
    const image_region& region) {
  if (!glimage) throw std::runtime_error("glimage is not initialized");
  update_gltexture_region(glimage.texture, img, region, glimage.texture.mipmap);
}
void update_glimage_region(opengl_image& glimage, const image<vec4b>& img,
    const image_region& region) {
  if (!glimage) throw std::runtime_error("glimage is not initialized");
  update_gltexture_region(glimage.texture, img, region, glimage.texture.mipmap);
}

// draw image
void draw_glimage(opengl_image& glimage, const draw_glimage_params& params) {
  check_glerror();
  set_glviewport(params.framebuffer);
  clear_glframebuffer(params.background);
  bind_glprogram(glimage.program);
  set_gluniform_texture(glimage.program, "txt", glimage.texture, 0);
  set_gluniform(glimage.program, "window_size",
      vec2f{(float)params.window.x, (float)params.window.y});
  set_gluniform(glimage.program, "image_size",
      vec2f{(float)glimage.texture.size.x, (float)glimage.texture.size.y});
  set_gluniform(glimage.program, "image_center", params.center);
  set_gluniform(glimage.program, "image_scale", params.scale);
  set_glvertexattrib(glimage.program, "texcoord", glimage.texcoord, zero2f);
  draw_gltriangles(glimage.element, 2);
  unbind_opengl_program();
  check_glerror();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize an OpenGL scene
void make_glscene(opengl_scene& glscene) {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

  static const char* vertex =
      R"(
        #version 330

        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;                    // shape transform
        uniform mat4 shape_xform_invtranspose;       // shape transform
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
            norm = (shape_xform_invtranspose * vec4(norm,0)).xyz;
            tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
        }
        )";

  static const char* fragment =
      R"(
        #version 330

        float pif = 3.14159265;

        uniform bool eyelight;         // eyelight shading
        uniform vec3 lamb;             // ambient light
        uniform int lnum;              // number of lights
        uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
        uniform vec3 lpos[16];         // light positions
        uniform vec3 lke[16];          // light intensities

        void evaluate_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
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

        vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
            vec3 n, vec3 wi, vec3 wo) {
            if(etype == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(rs*rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(etype == 1) {
                return ((1+dot(wo,wi))/2) * kd/pif;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pif;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pif)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pif;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pif)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        uniform int elem_type;
        uniform bool elem_faceted;
        uniform vec4 highlight;   // highlighted color

        uniform int mat_type;          // material type
        uniform vec3 mat_ke;           // material ke
        uniform vec3 mat_kd;           // material kd
        uniform vec3 mat_ks;           // material ks
        uniform float mat_rs;          // material rs
        uniform float mat_op;          // material op

        uniform bool mat_ke_txt_on;    // material ke texture on
        uniform sampler2D mat_ke_txt;  // material ke texture
        uniform bool mat_kd_txt_on;    // material kd texture on
        uniform sampler2D mat_kd_txt;  // material kd texture
        uniform bool mat_ks_txt_on;    // material ks texture on
        uniform sampler2D mat_ks_txt;  // material ks texture
        uniform bool mat_rs_txt_on;    // material rs texture on
        uniform sampler2D mat_rs_txt;  // material rs texture
        uniform bool mat_op_txt_on;    // material op texture on
        uniform sampler2D mat_op_txt;  // material op texture

        uniform bool mat_norm_txt_on;    // material norm texture on
        uniform sampler2D mat_norm_txt;  // material norm texture

        uniform bool mat_double_sided;   // double sided rendering

        uniform mat4 shape_xform;              // shape transform
        uniform mat4 shape_xform_invtranspose; // shape transform

        bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                           out vec3 kd, out vec3 ks, out float rs, out float op) {
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

            vec4 ke_txt = (mat_ke_txt_on) ? texture(mat_ke_txt,texcoord) : vec4(1,1,1,1);
            vec4 kd_txt = (mat_kd_txt_on) ? texture(mat_kd_txt,texcoord) : vec4(1,1,1,1);
            vec4 ks_txt = (mat_ks_txt_on) ? texture(mat_ks_txt,texcoord) : vec4(1,1,1,1);
            vec4 rs_txt = (mat_rs_txt_on) ? texture(mat_rs_txt,texcoord) : vec4(1,1,1,1);
            vec4 op_txt = (mat_op_txt_on) ? texture(mat_op_txt,texcoord) : vec4(1,1,1,1);

            // get material color from textures and adjust values
            if(mat_type == 1) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt.y;
                rs = rs*rs;
                op *= op_txt.x * kd_txt.w;
            } else if(mat_type == 2) {
                ke *= ke_txt.xyz;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(mat_type == 3) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                float gs = (1 - rs) * ks_txt.w;
                rs = 1 - gs;
                rs = rs*rs;
                op *= kd_txt.w;
            }

            return true;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!mat_norm_txt_on) return norm;
            vec3 tangu = normalize((shape_xform * vec4(normalize(tangsp.xyz),0)).xyz);
            vec3 tangv = normalize(cross(norm, tangu));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 texture = 2 * texture(mat_norm_txt,texcoord).xyz - 1;
            texture.y = -texture.y;
            return normalize( tangu * texture.x + tangv * texture.y + norm * texture.z );
        }

        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        uniform vec3 cam_pos;          // camera position
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        uniform float exposure; 
        uniform float gamma;

        out vec4 frag_color;      

        vec3 triangle_normal(vec3 pos) {
            vec3 fdx = dFdx(pos); 
            vec3 fdy = dFdy(pos); 
            return normalize((shape_xform * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
        }

        // main
        void main() {
            // view vector
            vec3 wo = normalize(cam_pos - pos);

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
            vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
            bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

            // exit if needed
            if(brdf_op < 0.005) discard;

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
                    c += pif * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lamb * brdf_kd;
                    // foreach light
                    for(int lid = 0; lid < lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        evaluate_light(lid, pos, cl, wi);
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

  // load program
  init_glprogram(glscene.program, vertex, fragment);
}

// Draw a shape
void draw_glinstance(opengl_scene& state, const opengl_instance& instance,
    const draw_glscene_params& params) {
  auto& shape    = state.shapes[instance.shape];
  auto& material = state.materials[instance.material];

  set_gluniform(state.program, "shape_xform", mat4f(instance.frame));
  set_gluniform(state.program, "shape_xform_invtranspose",
      transpose(mat4f(inverse(instance.frame, params.non_rigid_frames))));
  set_gluniform(state.program, "shape_normal_offset", 0.0f);
  set_gluniform(state.program, "highlight",
      instance.highlighted ? vec4f{1, 1, 0, 1} : zero4f);

  auto mtype = 2;
  if (material.gltf_textures) mtype = 3;
  set_gluniform(state.program, "mat_type", mtype);
  set_gluniform(state.program, "mat_ke", material.emission);
  set_gluniform(state.program, "mat_kd", material.diffuse);
  set_gluniform(state.program, "mat_ks", vec3f{material.metallic});
  set_gluniform(state.program, "mat_rs", material.roughness);
  set_gluniform(state.program, "mat_op", material.opacity);
  set_gluniform(state.program, "mat_double_sided", (int)params.double_sided);
  if (material.emission_map >= 0) {
    set_gluniform_texture(state.program, "mat_ke_txt", "mat_ke_txt_on",
        state.textures.at(material.emission_map), 0);
  } else {
    set_gluniform_texture(
        state.program, "mat_ke_txt", "mat_ke_txt_on", opengl_texture{}, 0);
  }
  if (material.diffuse_map >= 0) {
    set_gluniform_texture(state.program, "mat_kd_txt", "mat_kd_txt_on",
        state.textures.at(material.diffuse_map), 1);
  } else {
    set_gluniform_texture(
        state.program, "mat_kd_txt", "mat_kd_txt_on", opengl_texture{}, 1);
  }
  if (material.metallic_map >= 0) {
    set_gluniform_texture(state.program, "mat_ks_txt", "mat_ks_txt_on",
        state.textures.at(material.metallic_map), 2);
  } else {
    set_gluniform_texture(
        state.program, "mat_ks_txt", "mat_ks_txt_on", opengl_texture{}, 2);
  }
  if (material.roughness_map >= 0) {
    set_gluniform_texture(state.program, "mat_rs_txt", "mat_rs_txt_on",
        state.textures.at(material.roughness_map), 3);
  } else {
    set_gluniform_texture(
        state.program, "mat_rs_txt", "mat_rs_txt_on", opengl_texture{}, 3);
  }
  if (material.normal_map >= 0) {
    set_gluniform_texture(state.program, "mat_norm_txt", "mat_norm_txt_on",
        state.textures.at(material.normal_map), 5);
  } else {
    set_gluniform_texture(
        state.program, "mat_norm_txt", "mat_norm_txt_on", opengl_texture{}, 5);
  }

  set_gluniform(state.program, "elem_faceted", (int)!shape.normals);
  set_glvertexattrib(state.program, "vert_pos", shape.positions, zero3f);
  set_glvertexattrib(state.program, "vert_norm", shape.normals, zero3f);
  set_glvertexattrib(state.program, "vert_texcoord", shape.texcoords, zero2f);
  set_glvertexattrib(
      state.program, "vert_color", shape.colors, vec4f{1, 1, 1, 1});
  set_glvertexattrib(
      state.program, "vert_tangsp", shape.tangentsps, vec4f{0, 0, 1, 1});

  if (shape.points) {
    set_gluniform(state.program, "elem_type", 1);
    draw_glpoints(shape.points, shape.points.num);
  }
  if (shape.lines) {
    set_gluniform(state.program, "elem_type", 2);
    draw_gllines(shape.lines, shape.lines.num);
  }
  if (shape.triangles) {
    set_gluniform(state.program, "elem_type", 3);
    draw_gltriangles(shape.triangles, shape.triangles.num);
  }
  if (shape.quads) {
    set_gluniform(state.program, "elem_type", 3);
    draw_gltriangles(shape.quads, shape.quads.num);
  }

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_glculling(false);
        check_glerror();
        set_gluniform(state.program, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state.program, "ke"), 0, 0, 0);
        set_gluniform(state.program, "op"), material.op);
        set_gluniform(state.program, "shp_normal_offset"), 0.01f);
        check_glerror();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_glerror();
    }
#endif
  if (params.edges) throw std::runtime_error("edges are momentarily disabled");

  // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
}

// Display a scene
void draw_glscene(opengl_scene& state, const vec4i& viewport,
    const draw_glscene_params& params) {
  auto& glcamera    = state.cameras.at(params.camera);
  auto  camera_view = mat4f(inverse(glcamera.frame));
  auto  camera_proj = perspective_mat(glcamera.yfov,
      (float)viewport.z / (float)viewport.w, params.near, params.far);

  clear_glframebuffer(params.background);
  set_glviewport(viewport);

  bind_glprogram(state.program);
  set_gluniform(state.program, "cam_pos", glcamera.frame.o);
  set_gluniform(state.program, "cam_xform_inv", camera_view);
  set_gluniform(state.program, "cam_proj", camera_proj);
  set_gluniform(state.program, "eyelight", (int)params.eyelight);
  set_gluniform(state.program, "exposure", params.exposure);
  set_gluniform(state.program, "gamma", params.gamma);

  if (!params.eyelight) {
    set_gluniform(state.program, "lamb", zero3f);
    set_gluniform(state.program, "lnum", (int)state.lights.size());
    for (auto i = 0; i < state.lights.size(); i++) {
      auto is = std::to_string(i);
      set_gluniform(state.program, ("lpos[" + is + "]").c_str(),
          state.lights[i].position);
      set_gluniform(
          state.program, ("lke[" + is + "]").c_str(), state.lights[i].emission);
      set_gluniform(state.program, ("ltype[" + is + "]").c_str(),
          (int)state.lights[i].type);
    }
  }

  if (params.wireframe) set_glwireframe(true);
  for (auto& instance : state.instances) {
    draw_glinstance(state, instance, params);
  }

  unbind_opengl_program();
  if (params.wireframe) set_glwireframe(false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

opengl_program::opengl_program(opengl_program&& other) {
  operator=(std::forward<opengl_program>(other));
}
opengl_program& opengl_program::operator=(opengl_program&& other) {
  std::swap(program_id, other.program_id);
  std::swap(vertex_shader_id, other.vertex_shader_id);
  std::swap(fragment_shader_id, other.fragment_shader_id);
  std::swap(vertex_array_object_id, other.vertex_array_object_id);
  return *this;
}
opengl_program::~opengl_program() { delete_glprogram(*this); }

void init_glprogram(
    opengl_program& program, const char* vertex, const char* fragment) {
  assert(glGetError() == GL_NO_ERROR);
  glGenVertexArrays(1, &program.vertex_array_object_id);
  glBindVertexArray(program.vertex_array_object_id);
  assert(glGetError() == GL_NO_ERROR);

  int  errflags;
  char errbuf[10000];

  // create vertex
  assert(glGetError() == GL_NO_ERROR);
  program.vertex_shader_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(program.vertex_shader_id, 1, &vertex, NULL);
  glCompileShader(program.vertex_shader_id);
  glGetShaderiv(program.vertex_shader_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program.vertex_shader_id, 10000, 0, errbuf);
    throw std::runtime_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create fragment
  assert(glGetError() == GL_NO_ERROR);
  program.fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(program.fragment_shader_id, 1, &fragment, NULL);
  glCompileShader(program.fragment_shader_id);
  glGetShaderiv(program.fragment_shader_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(program.fragment_shader_id, 10000, 0, errbuf);
    throw std::runtime_error("shader not compiled with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);

  // create program
  assert(glGetError() == GL_NO_ERROR);
  program.program_id = glCreateProgram();
  glAttachShader(program.program_id, program.vertex_shader_id);
  glAttachShader(program.program_id, program.fragment_shader_id);
  glLinkProgram(program.program_id);
  glValidateProgram(program.program_id);
  glGetProgramiv(program.program_id, GL_LINK_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
    throw std::runtime_error("program not linked with error\n"s + errbuf);
  }
  glGetProgramiv(program.program_id, GL_VALIDATE_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program.program_id, 10000, 0, errbuf);
    throw std::runtime_error("program not linked with error\n"s + errbuf);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void delete_glprogram(opengl_program& program) {
  if (!program) return;
  glDeleteProgram(program.program_id);
  glDeleteShader(program.vertex_shader_id);
  glDeleteShader(program.fragment_shader_id);
  program.program_id         = 0;
  program.vertex_shader_id   = 0;
  program.fragment_shader_id = 0;
}

opengl_texture::opengl_texture(opengl_texture&& other) {
  operator=(std::forward<opengl_texture>(other));
}
opengl_texture& opengl_texture::operator=(opengl_texture&& other) {
  std::swap(texture_id, other.texture_id);
  std::swap(size, other.size);
  return *this;
}
opengl_texture::~opengl_texture() { delete_gltexture(*this); }

void init_gltexture(opengl_texture& texture, const vec2i& size, bool as_float,
    bool as_srgb, bool linear, bool mipmap) {
  if (texture) delete_gltexture(texture);
  assert(glGetError() == GL_NO_ERROR);
  glGenTextures(1, &texture.texture_id);
  texture.size     = size;
  texture.mipmap   = mipmap;
  texture.is_srgb  = as_srgb;
  texture.is_float = as_float;
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  if (as_float) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  } else if (as_srgb) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, size.x, size.y, 0, GL_RGBA,
        GL_UNSIGNED_BYTE, nullptr);
  } else {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  }
  if (mipmap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture(
    opengl_texture& texture, const image<vec4f>& img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
      GL_FLOAT, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture_region(opengl_texture& texture, const image<vec4f>& img,
    const image_region& region, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  auto clipped = image<vec4f>{};
  get_region(clipped, img, region);
  glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y, region.size().x,
      region.size().y, GL_RGBA, GL_FLOAT, clipped.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture(
    opengl_texture& texture, const image<vec4b>& img, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.size().x, img.size().y, GL_RGBA,
      GL_UNSIGNED_BYTE, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void update_gltexture_region(opengl_texture& texture, const image<vec4b>& img,
    const image_region& region, bool mipmap) {
  assert(glGetError() == GL_NO_ERROR);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  auto clipped = image<vec4b>{};
  get_region(clipped, img, region);
  glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y, region.size().x,
      region.size().y, GL_RGBA, GL_UNSIGNED_BYTE, clipped.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  assert(glGetError() == GL_NO_ERROR);
}

void delete_gltexture(opengl_texture& texture) {
  if (!texture) return;
  glDeleteTextures(1, &texture.texture_id);
  texture.texture_id = 0;
  texture.size       = zero2i;
}

opengl_arraybuffer::opengl_arraybuffer(opengl_arraybuffer&& other) {
  operator=(std::forward<opengl_arraybuffer>(other));
}
opengl_arraybuffer& opengl_arraybuffer::operator=(opengl_arraybuffer&& other) {
  std::swap(buffer_id, other.buffer_id);
  std::swap(num, other.num);
  std::swap(elem_size, elem_size);
  return *this;
}
opengl_arraybuffer::~opengl_arraybuffer() { delete_glarraybuffer(*this); }

template <typename T>
void init_glarray_buffer_impl(
    opengl_arraybuffer& buffer, const vector<T>& array, bool dynamic) {
  buffer           = opengl_arraybuffer{};
  buffer.num       = size(array);
  buffer.elem_size = sizeof(T);
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer.buffer_id);
  glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
  glBufferData(GL_ARRAY_BUFFER, size(array) * sizeof(T), array.data(),
      (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<float>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec2f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec3f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}
void init_glarraybuffer(
    opengl_arraybuffer& buffer, const vector<vec4f>& data, bool dynamic) {
  init_glarray_buffer_impl(buffer, data, dynamic);
}

void delete_glarraybuffer(opengl_arraybuffer& buffer) {
  if (!buffer) return;
  glDeleteBuffers(1, &buffer.buffer_id);
  buffer.buffer_id = 0;
  buffer.elem_size = 0;
  buffer.num       = 0;
}

opengl_elementbuffer::opengl_elementbuffer(opengl_elementbuffer&& other) {
  operator=(std::forward<opengl_elementbuffer>(other));
}
opengl_elementbuffer& opengl_elementbuffer::operator=(
    opengl_elementbuffer&& other) {
  std::swap(buffer_id, other.buffer_id);
  std::swap(num, other.num);
  std::swap(elem_size, elem_size);
  return *this;
}
opengl_elementbuffer::~opengl_elementbuffer() { delete_glelementbuffer(*this); }

template <typename T>
void init_glelementbuffer_impl(
    opengl_elementbuffer& buffer, const vector<T>& array, bool dynamic) {
  buffer           = opengl_elementbuffer{};
  buffer.num       = size(array);
  buffer.elem_size = sizeof(T);
  assert(glGetError() == GL_NO_ERROR);
  glGenBuffers(1, &buffer.buffer_id);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, size(array) * sizeof(T), array.data(),
      (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  assert(glGetError() == GL_NO_ERROR);
}

void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<int>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}
void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<vec2i>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}
void init_glelementbuffer(
    opengl_elementbuffer& buffer, const vector<vec3i>& data, bool dynamic) {
  init_glelementbuffer_impl(buffer, data, dynamic);
}

void delete_glelementbuffer(opengl_elementbuffer& buffer) {
  if (!buffer) return;
  glDeleteBuffers(1, &buffer.buffer_id);
  buffer.buffer_id = 0;
  buffer.elem_size = 0;
  buffer.num       = 0;
}

void bind_glprogram(opengl_program& program) {
  glUseProgram(program.program_id);
}
void unbind_opengl_program() { glUseProgram(0); }

int get_gluniform_location(const opengl_program& program, const char* name) {
  return glGetUniformLocation(program.program_id, name);
}

void set_gluniform(int locatiom, int value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform1i(locatiom, value);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec2i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform2i(locatiom, value.x, value.y);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec3i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform3i(locatiom, value.x, value.y, value.z);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec4i& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform4i(locatiom, value.x, value.y, value.z, value.w);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, float value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform1f(locatiom, value);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec2f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform2f(locatiom, value.x, value.y);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform3f(locatiom, value.x, value.y, value.z);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const vec4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniform4f(locatiom, value.x, value.y, value.z, value.w);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const mat4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniformMatrix4fv(locatiom, 1, false, &value.x.x);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform(int locatiom, const frame3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  glUniformMatrix4x3fv(locatiom, 1, false, &value.x.x);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(
    int locatiom, const opengl_texture& texture, int unit) {
  assert(glGetError() == GL_NO_ERROR);
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture.texture_id);
  glUniform1i(locatiom, unit);
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(opengl_program& program, const char* name,
    const opengl_texture& texture, int unit) {
  set_gluniform_texture(get_gluniform_location(program, name), texture, unit);
}

void set_gluniform_texture(
    int locatiom, int locatiom_on, const opengl_texture& texture, int unit) {
  assert(glGetError() == GL_NO_ERROR);
  if (texture.texture_id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.texture_id);
    glUniform1i(locatiom, unit);
    glUniform1i(locatiom_on, 1);
  } else {
    glUniform1i(locatiom_on, 0);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_gluniform_texture(opengl_program& program, const char* name,
    const char* name_on, const opengl_texture& texture, int unit) {
  set_gluniform_texture(get_gluniform_location(program, name),
      get_gluniform_location(program, name_on), texture, unit);
}

int get_glvertexattrib_location(
    const opengl_program& program, const char* name) {
  return glGetAttribLocation(program.program_id, name);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, float value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 1, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib1f(locatiom, value);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec2f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 2, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib2f(locatiom, value.x, value.y);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec3f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 3, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib3f(locatiom, value.x, value.y, value.z);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec4f& value) {
  assert(glGetError() == GL_NO_ERROR);
  if (buffer.buffer_id) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer.buffer_id);
    glEnableVertexAttribArray(locatiom);
    glVertexAttribPointer(locatiom, 4, GL_FLOAT, false, 0, nullptr);
  } else {
    glVertexAttrib4f(locatiom, value.x, value.y, value.z, value.w);
  }
  assert(glGetError() == GL_NO_ERROR);
}

void draw_glpoints(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_POINTS, num, GL_UNSIGNED_INT, nullptr);
}

void draw_gllines(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_LINES, num * 2, GL_UNSIGNED_INT, nullptr);
}

void draw_gltriangles(const opengl_elementbuffer& buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.buffer_id);
  glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

void draw_glimage(const opengl_texture& texture, int win_width, int win_height,
    const vec2f& image_center, float image_scale) {
  static opengl_program       gl_prog      = {};
  static opengl_arraybuffer   gl_texcoord  = {};
  static opengl_elementbuffer gl_triangles = {};

  // initialization
  if (!gl_prog) {
    auto vert = R"(
            #version 330
            in vec2 texcoord;
            out vec2 frag_texcoord;
            uniform vec2 window_size, image_size;
            uniform vec2 image_center;
            uniform float image_scale;
            void main() {
                vec2 pos = (texcoord - vec2(0.5,0.5)) * image_size * image_scale + image_center;
                gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
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
    init_glprogram(gl_prog, vert, frag);
    init_glarraybuffer(
        gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    init_glelementbuffer(
        gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
  }

  // draw
  check_glerror();
  bind_glprogram(gl_prog);
  set_gluniform_texture(gl_prog, "txt", texture, 0);
  set_gluniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_gluniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_gluniform(gl_prog, "image_center", image_center);
  set_gluniform(gl_prog, "image_scale", image_scale);
  set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
  draw_gltriangles(gl_triangles, 2);
  unbind_opengl_program();
  check_glerror();
}

void draw_glimage_background(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size) {
  static opengl_program       gl_prog      = {};
  static opengl_arraybuffer   gl_texcoord  = {};
  static opengl_elementbuffer gl_triangles = {};

  // initialization
  if (!gl_prog) {
    auto vert = R"(
            #version 330
            in vec2 texcoord;
            out vec2 frag_texcoord;
            uniform vec2 window_size, image_size, border_size;
            uniform vec2 image_center;
            uniform float image_scale;
            void main() {
                vec2 pos = (texcoord - vec2(0.5,0.5)) * (image_size + border_size*2) * image_scale + image_center;
                gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
                frag_texcoord = texcoord;
            }
        )";
    auto frag = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform vec2 image_size, border_size;
            uniform float image_scale;
            void main() {
                ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
                ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
                ivec2 tile = tilecoord / 16;
                if(imcoord.x <= 0 || imcoord.y <= 0 || 
                    imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
                else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
                else frag_color = vec4(0.3,0.3,0.3,1);
            }
        )";
    init_glprogram(gl_prog, vert, frag);
    init_glarraybuffer(
        gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    init_glelementbuffer(
        gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
  }

  // draw
  bind_glprogram(gl_prog);
  set_gluniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_gluniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_gluniform(
      gl_prog, "border_size", vec2f{(float)border_size, (float)border_size});
  set_gluniform(gl_prog, "image_center", image_center);
  set_gluniform(gl_prog, "image_scale", image_scale);
  set_glvertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
  draw_gltriangles(gl_triangles, 2);
  unbind_opengl_program();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

void _glfw_refresh_callback(GLFWwindow* glfw) {
  auto& win = *(const opengl_window*)glfwGetWindowUserPointer(glfw);
  if (win.refresh_cb) win.refresh_cb(win);
}

void _glfw_drop_callback(GLFWwindow* glfw, int num, const char** paths) {
  auto& win = *(const opengl_window*)glfwGetWindowUserPointer(glfw);
  if (win.drop_cb) {
    auto pathv = vector<string>();
    for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
    win.drop_cb(win, pathv);
  }
}

void init_glwindow(opengl_window& win, const vec2i& size, const string& title,
    void* user_pointer, refresh_glcallback refresh_cb) {
  // init glfw
  if (!glfwInit())
    throw std::runtime_error("cannot initialize windowing system");
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // create window
  win     = opengl_window();
  win.win = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
  if (!win.win) throw std::runtime_error("cannot initialize windowing system");
  glfwMakeContextCurrent(win.win);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowRefreshCallback(win.win, _glfw_refresh_callback);
  glfwSetWindowUserPointer(win.win, &win);
  win.user_ptr   = user_pointer;
  win.refresh_cb = refresh_cb;

  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error("cannot initialize OpenGL extensions");
}

void delete_glwindow(opengl_window& win) {
  glfwDestroyWindow(win.win);
  glfwTerminate();
  win.win = nullptr;
}

void* get_gluser_pointer(const opengl_window& win) { return win.user_ptr; }

void set_drop_glcallback(opengl_window& win, drop_glcallback drop_cb) {
  win.drop_cb = drop_cb;
  glfwSetDropCallback(win.win, _glfw_drop_callback);
}

vec2i get_glframebuffer_size(const opengl_window& win, bool ignore_widgets) {
  auto size = zero2i;
  glfwGetFramebufferSize(win.win, &size.x, &size.y);
  if (ignore_widgets && win.widgets_width) {
    auto win_size = zero2i;
    glfwGetWindowSize(win.win, &win_size.x, &win_size.y);
    size.x -= (int)(win.widgets_width * (float)size.x / (float)win_size.x);
  }
  return size;
}

vec4i get_glframebuffer_viewport(
    const opengl_window& win, bool ignore_widgets) {
  auto viewport = zero4i;
  glfwGetFramebufferSize(win.win, &viewport.z, &viewport.w);
  if (ignore_widgets && win.widgets_width) {
    auto win_size = zero2i;
    glfwGetWindowSize(win.win, &win_size.x, &win_size.y);
    auto offset = (int)(win.widgets_width * (float)viewport.z / win_size.x);
    viewport.z -= offset;
    if (win.widgets_left) viewport.x += offset;
  }
  return viewport;
}

vec2i get_glwindow_size(const opengl_window& win, bool ignore_widgets) {
  auto size = zero2i;
  glfwGetWindowSize(win.win, &size.x, &size.y);
  if (ignore_widgets && win.widgets_width) size.x -= win.widgets_width;
  return size;
}

bool should_glwindow_close(const opengl_window& win) {
  return glfwWindowShouldClose(win.win);
}
void set_glwindow_close(const opengl_window& win, bool close) {
  glfwSetWindowShouldClose(win.win, close ? GLFW_TRUE : GLFW_FALSE);
}

vec2f get_glmouse_pos(const opengl_window& win, bool ignore_widgets) {
  double mouse_posx, mouse_posy;
  glfwGetCursorPos(win.win, &mouse_posx, &mouse_posy);
  auto pos = vec2f{(float)mouse_posx, (float)mouse_posy};
  if (ignore_widgets && win.widgets_width && win.widgets_left) {
    pos.x -= win.widgets_width;
  }
  return pos;
}

bool get_glmouse_left(const opengl_window& win) {
  return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
}
bool get_glmouse_right(const opengl_window& win) {
  return glfwGetMouseButton(win.win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
}

bool get_glalt_key(const opengl_window& win) {
  return glfwGetKey(win.win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
         glfwGetKey(win.win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

bool get_glshift_key(const opengl_window& win) {
  return glfwGetKey(win.win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
         glfwGetKey(win.win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

void process_glevents(const opengl_window& win, bool wait) {
  if (wait)
    glfwWaitEvents();
  else
    glfwPollEvents();
}

void swap_glbuffers(const opengl_window& win) { glfwSwapBuffers(win.win); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

void init_glwidgets(opengl_window& win, int width, bool left) {
  // init widgets
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(win.win, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();
  win.widgets_width = width;
  win.widgets_left  = left;
}

bool get_glwidgets_active(const opengl_window& win) {
  auto io = &ImGui::GetIO();
  return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

void begin_glwidgets(const opengl_window& win) {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  auto win_size = get_glwindow_size(win, false);
  if (win.widgets_left) {
    ImGui::SetNextWindowPos({0, 0});
    ImGui::SetNextWindowSize({(float)win.widgets_width, (float)win_size.y});
  } else {
    ImGui::SetNextWindowPos({(float)(win_size.x - win.widgets_width), 0});
    ImGui::SetNextWindowSize({(float)win.widgets_width, (float)win_size.y});
  }
  ImGui::SetNextWindowCollapsed(false);
  ImGui::SetNextWindowBgAlpha(1);
}

void end_glwidgets(const opengl_window& win) {
  ImGui::End();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool begin_glwidgets_window(const opengl_window& win, const char* title) {
  return ImGui::Begin(title, nullptr,
      // ImGuiWindowFlags_NoTitleBar |
      ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
          ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings);
}

bool begin_glheader(const opengl_window& win, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_glheader(const opengl_window& win) { ImGui::PopID(); }

void open_glmodal(const opengl_window& win, const char* lbl) {
  ImGui::OpenPopup(lbl);
}
void clear_glmodal(const opengl_window& win) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(const opengl_window& win, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(const opengl_window& win) { ImGui::EndPopup(); }
bool is_glmodal_open(const opengl_window& win, const char* lbl) {
  return ImGui::IsPopupOpen(lbl);
}

bool draw_glmessage(
    const opengl_window& win, const char* lbl, const string& message) {
  if (ImGui::BeginPopupModal(lbl)) {
    auto open = true;
    ImGui::Text("%s", message.c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      open = false;
    }
    ImGui::EndPopup();
    return open;
  } else {
    return false;
  }
}

std::deque<string> _message_queue = {};
std::mutex         _message_mutex;
void               push_glmessage(const string& message) {
  std::lock_guard lock(_message_mutex);
  _message_queue.push_back(message);
}
void push_glmessage(const opengl_window& win, const string& message) {
  std::lock_guard lock(_message_mutex);
  _message_queue.push_back(message);
}
bool draw_glmessages(const opengl_window& win) {
  std::lock_guard lock(_message_mutex);
  if (_message_queue.empty()) return false;
  if (!is_glmodal_open(win, "<message>")) {
    open_glmodal(win, "<message>");
    return true;
  } else if (ImGui::BeginPopupModal("<message>")) {
    ImGui::Text("%s", _message_queue.front().c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      _message_queue.pop_front();
    }
    ImGui::EndPopup();
    return true;
  } else {
    return false;
  }
}

struct filedialog_state {
  string                     dirname       = "";
  string                     filename      = "";
  vector<pair<string, bool>> entries       = {};
  bool                       save          = false;
  bool                       remove_hidden = true;
  string                     filter        = "";
  vector<string>             extensions    = {};

  filedialog_state() {}
  filedialog_state(const string& dirname, const string& filename, bool save,
      const string& filter) {
    this->save = save;
    set_filter(filter);
    set_dirname(dirname);
    set_filename(filename);
  }
  void set_dirname(const string& name) {
    dirname = name;
    dirname = normalize_path(dirname);
    if (dirname == "") dirname = "./";
    if (dirname.back() != '/') dirname += '/';
    refresh();
  }
  void set_filename(const string& name) {
    filename = name;
    check_filename();
  }
  void set_filter(const string& flt) {
    auto globs = vector<string>{""};
    for (auto i = 0; i < flt.size(); i++) {
      if (flt[i] == ';') {
        globs.push_back("");
      } else {
        globs.back() += flt[i];
      }
    }
    filter = "";
    extensions.clear();
    for (auto pattern : globs) {
      if (pattern == "") continue;
      auto ext = get_extension(pattern);
      if (ext != "") {
        extensions.push_back(ext);
        filter += (filter == "") ? ("*." + ext) : (";*." + ext);
      }
    }
  }
  void check_filename() {
    if (filename.empty()) return;
    auto ext = get_extension(filename);
    if (std::find(extensions.begin(), extensions.end(), ext) ==
        extensions.end()) {
      filename = "";
      return;
    }
    if (!save && !exists_file(dirname + filename)) {
      filename = "";
      return;
    }
  }
  void select_entry(int idx) {
    if (entries[idx].second) {
      set_dirname(dirname + entries[idx].first);
    } else {
      set_filename(entries[idx].first);
    }
  }

  void refresh() {
    entries.clear();
    cf_dir_t dir;
    cf_dir_open(&dir, dirname.c_str());
    while (dir.has_next) {
      cf_file_t file;
      cf_read_file(&dir, &file);
      cf_dir_next(&dir);
      if (remove_hidden && file.name[0] == '.') continue;
      if (file.is_dir) {
        entries.push_back({file.name + "/"s, true});
      } else {
        entries.push_back({file.name, false});
      }
    }
    cf_dir_close(&dir);
    std::sort(entries.begin(), entries.end(), [](auto& a, auto& b) {
      if (a.second == b.second) return a.first < b.first;
      return a.second;
    });
  }

  string get_path() const { return dirname + filename; }
  bool   exists_file(const string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
  }
};
bool draw_glfiledialog(const opengl_window& win, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter) {
  static auto states = hash_map<string, filedialog_state>{};
  ImGui::SetNextWindowSize({500, 300}, ImGuiCond_FirstUseEver);
  if (ImGui::BeginPopupModal(lbl)) {
    if (states.find(lbl) == states.end()) {
      states[lbl] = filedialog_state{dirname, filename, save, filter};
    }
    auto& state = states.at(lbl);
    char  dir_buffer[1024];
    strcpy(dir_buffer, state.dirname.c_str());
    if (ImGui::InputText("dir", dir_buffer, sizeof(dir_buffer))) {
      state.set_dirname(dir_buffer);
    }
    auto current_item = -1;
    if (ImGui::ListBox(
            "entries", &current_item,
            [](void* data, int idx, const char** out_text) -> bool {
              auto& state = *(filedialog_state*)data;
              *out_text   = state.entries[idx].first.c_str();
              return true;
            },
            &state, (int)state.entries.size())) {
      state.select_entry(current_item);
    }
    char file_buffer[1024];
    strcpy(file_buffer, state.filename.c_str());
    if (ImGui::InputText("file", file_buffer, sizeof(file_buffer))) {
      state.set_filename(file_buffer);
    }
    char filter_buffer[1024];
    strcpy(filter_buffer, state.filter.c_str());
    if (ImGui::InputText("filter", filter_buffer, sizeof(filter_buffer))) {
      state.set_filter(filter_buffer);
    }
    auto ok = false, exit = false;
    if (ImGui::Button("Ok")) {
      path = state.dirname + state.filename;
      ok   = true;
      exit = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      exit = true;
    }
    if (exit) {
      ImGui::CloseCurrentPopup();
      states.erase(lbl);
    }
    ImGui::EndPopup();
    return ok;
  } else {
    return false;
  }
}
bool draw_glfiledialog_button(const opengl_window& win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter) {
  if (is_glmodal_open(win, lbl)) {
    return draw_glfiledialog(win, lbl, path, save, dirname, filename, filter);
  } else {
    if (draw_glbutton(win, button_lbl, button_active)) {
      open_glmodal(win, lbl);
    }
    return false;
  }
}

bool draw_glbutton(const opengl_window& win, const char* lbl) {
  return ImGui::Button(lbl);
}
bool draw_glbutton(const opengl_window& win, const char* lbl, bool enabled) {
  if (enabled) {
    return ImGui::Button(lbl);
  } else {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    auto ok = ImGui::Button(lbl);
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
    return ok;
  }
}

void draw_gllabel(
    const opengl_window& win, const char* lbl, const string& texture) {
  ImGui::LabelText(lbl, "%s", texture.c_str());
}

void draw_glseparator(const opengl_window& win) { ImGui::Separator(); }

void continue_glline(const opengl_window& win) { ImGui::SameLine(); }

bool draw_gltextinput(
    const opengl_window& win, const char* lbl, string& value) {
  char buffer[4096];
  auto num = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer, sizeof(buffer));
  if (edited) value = buffer;
  return edited;
}

bool draw_glslider(const opengl_window& win, const char* lbl, float& value,
    float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec2f& value,
    float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec3f& value,
    float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_glslider(const opengl_window& win, const char* lbl, vec4f& value,
    float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_glslider(
    const opengl_window& win, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_gldragger(const opengl_window& win, const char* lbl, float& value,
    float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_gldragger(const opengl_window& win, const char* lbl, int& value,
    float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2i& value,
    float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3i& value,
    float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4i& value,
    float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_glcheckbox(const opengl_window& win, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}

bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(xyz(color));
  if (scale > 1) {
    xyz(color) /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    xyz(value) = xyz(color) * exp2(exposure);
    value.w    = color.w;
    return true;
  } else {
    return false;
  }
}

bool draw_glcombobox(const opengl_window& win, const char* lbl, int& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, labels[value].c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(const opengl_window& win, const char* lbl, string& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(const opengl_window& win, const char* lbl, int& idx,
    int num, const std::function<const char*(int)>& labels, bool include_null) {
  if (num <= 0) idx = -1;
  if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx) : "<none>")) return false;
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

void draw_glhistogram(
    const opengl_window& win, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec4f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " w"s).c_str(), (const float*)values.data() + 3,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
}

// https://github.com/ocornut/imgui/issues/300
struct ImGuiAppLog {
  ImGuiTextBuffer Buf;
  ImGuiTextFilter Filter;
  ImVector<int>   LineOffsets;  // Index to lines offset
  bool            ScrollToBottom;

  void Clear() {
    Buf.clear();
    LineOffsets.clear();
  }

  void AddLog(const char* msg, const char* lbl) {
    int old_size = Buf.size();
    Buf.appendf("[%s] %s\n", lbl, msg);
    for (int new_size = Buf.size(); old_size < new_size; old_size++)
      if (Buf[old_size] == '\n') LineOffsets.push_back(old_size);
    ScrollToBottom = true;
  }

  void Draw() {
    if (ImGui::Button("Clear")) Clear();
    ImGui::SameLine();
    bool copy = ImGui::Button("Copy");
    ImGui::SameLine();
    Filter.Draw("Filter", -100.0f);
    ImGui::Separator();
    ImGui::BeginChild("scrolling");
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
    if (copy) ImGui::LogToClipboard();

    if (Filter.IsActive()) {
      const char* buf_begin = Buf.begin();
      const char* line      = buf_begin;
      for (int line_no = 0; line != NULL; line_no++) {
        const char* line_end = (line_no < LineOffsets.Size)
                                   ? buf_begin + LineOffsets[line_no]
                                   : NULL;
        if (Filter.PassFilter(line, line_end))
          ImGui::TextUnformatted(line, line_end);
        line = line_end && line_end[1] ? line_end + 1 : NULL;
      }
    } else {
      ImGui::TextUnformatted(Buf.begin());
    }

    if (ScrollToBottom) ImGui::SetScrollHere(1.0f);
    ScrollToBottom = false;
    ImGui::PopStyleVar();
    ImGui::EndChild();
  }
  void Draw(const char* title, bool* p_opened = NULL) {
    ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(title, p_opened);
    Draw();
    ImGui::End();
  }
};

std::mutex  _log_mutex;
ImGuiAppLog _log_widget;
void        log_glinfo(const opengl_window& win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_glerror(const opengl_window& win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_gllogs(const opengl_window& win) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void draw_gllog(const opengl_window& win) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

}  // namespace yocto
