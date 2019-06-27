//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
#include "ysceneui.h"

#include <atomic>
#include <future>
#include <thread>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;
#ifdef _WIN32
#undef near
#undef far
#endif

#include "ext/CLI11.hpp"

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};

struct drawgl_shape {
  opengl_arraybuffer   positions_buffer     = {};
  opengl_arraybuffer   normals_buffer       = {};
  opengl_arraybuffer   texcoords_buffer     = {};
  opengl_arraybuffer   colors_buffer        = {};
  opengl_arraybuffer   tangentspaces_buffer = {};
  opengl_elementbuffer points_buffer        = {};
  opengl_elementbuffer lines_buffer         = {};
  opengl_elementbuffer triangles_buffer     = {};
  opengl_elementbuffer quads_buffer         = {};
};

struct drawgl_state {
  opengl_program         program  = {};
  vector<drawgl_shape>   shapes   = {};
  vector<opengl_texture> textures = {};
};

struct drawgl_lights {
  vector<vec3f> positions = {};
  vector<vec3f> emission  = {};
  vector<int>   types     = {};

  bool empty() const { return positions.empty(); }
};

void init_drawgl_lights(drawgl_lights& lights, const yocto_scene& scene) {
  lights = {};
  for (auto& instance : scene.instances) {
    if (instance.shape < 0) continue;
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    if (material.emission != zero3f) continue;
    if (lights.positions.size() >= 16) break;
    auto bbox = compute_bounds(shape);
    auto pos  = (bbox.max + bbox.min) / 2;
    auto area = 0.0f;
    if (!shape.triangles.empty()) {
      for (auto t : shape.triangles)
        area += triangle_area(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
    } else if (!shape.quads.empty()) {
      for (auto q : shape.quads)
        area += quad_area(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
    } else if (!shape.lines.empty()) {
      for (auto l : shape.lines)
        area += line_length(shape.positions[l.x], shape.positions[l.y]);
    } else {
      area += shape.positions.size();
    }
    auto ke = material.emission * area;
    lights.positions.push_back(transform_point(instance.frame, pos));
    lights.emission.push_back(ke);
    lights.types.push_back(0);
  }
}

// Draw options
struct drawgl_params {
  int   camera           = 0;
  int   resolution       = 1280;
  bool  wireframe        = false;
  bool  edges            = false;
  float edge_offset      = 0.01f;
  bool  eyelight         = false;
  float exposure         = 0;
  float gamma            = 2.2f;
  vec3f ambient          = {0, 0, 0};
  bool  double_sided     = true;
  bool  non_rigid_frames = true;
  float near             = 0.01f;
  float far              = 10000.0f;
};

// Equality operators
inline bool operator==(const drawgl_params& a, const drawgl_params& b) {
  return memcmp(&a, &b, sizeof(a)) == 0;
}
inline bool operator!=(const drawgl_params& a, const drawgl_params& b) {
  return memcmp(&a, &b, sizeof(a)) != 0;
}

// Application task
enum struct app_task_type {
  none,
  load_scene,
  load_element,
  apply_edit,
  save_image,
  save_scene,
  close_scene
};

struct app_task {
  app_task_type     type;
  std::future<void> result;
  std::atomic<bool> stop;
  std::atomic<int>  current;
  app_edit          edit;

  app_task(app_task_type type, const app_edit& edit = {})
      : type{type}, result{}, stop{false}, current{-1}, edit{edit} {}
  ~app_task() {
    stop = true;
    if (result.valid()) {
      try {
        result.get();
      } catch (...) {
      }
    }
  }
};

// Application state
struct app_scene {
  // loading parameters
  string filename  = "scene.json";
  string imagename = "out.png";
  string outname   = "scene.json";
  string name      = "";

  // options
  load_params   load_prms   = {};
  save_params   save_prms   = {};
  drawgl_params drawgl_prms = {};

  // scene
  yocto_scene scene = {};

  // rendering state
  drawgl_state  state  = {};
  drawgl_lights lights = {};

  // view image
  bool   navigation_fps = false;
  float  time           = 0;
  string anim_group     = "";
  vec2f  time_range     = zero2f;
  bool   animate        = false;

  // tasks
  bool                 load_done = false;
  std::deque<app_task> task_queue;
  app_selection        selection = {typeid(void), -1};
};

// Application state
struct app_state {
  // data
  std::deque<app_scene> scenes;
  int                   selected = -1;
  std::deque<string>    errors;

  // default options
  load_params   load_prms   = {};
  save_params   save_prms   = {};
  drawgl_params drawgl_prms = {};
};

void add_new_scene(app_state& app, const string& filename) {
  auto& scn       = app.scenes.emplace_back();
  scn.filename    = filename;
  scn.imagename   = fs::path(filename).replace_extension(".png");
  scn.outname     = fs::path(filename).replace_extension(".edited.yaml");
  scn.name        = fs::path(scn.filename).filename();
  scn.load_prms   = app.load_prms;
  scn.save_prms   = app.save_prms;
  scn.drawgl_prms = app.drawgl_prms;
  scn.task_queue.emplace_back(app_task_type::load_scene);
  app.selected = (int)app.scenes.size() - 1;
}

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

// Draw a shape
void draw_glinstance(drawgl_state& state, const yocto_scene& scene,
    const yocto_instance& instance, bool highlighted,
    const drawgl_params& options) {
  auto& shape    = scene.shapes[instance.shape];
  auto& vbos     = state.shapes.at(instance.shape);
  auto& material = scene.materials[instance.material];

  set_gluniform(state.program, "shape_xform", mat4f(instance.frame));
  set_gluniform(state.program, "shape_xform_invtranspose",
      transpose(mat4f(inverse(instance.frame, options.non_rigid_frames))));
  set_gluniform(state.program, "shape_normal_offset", 0.0f);
  set_gluniform(
      state.program, "highlight", (highlighted) ? vec4f{1, 1, 0, 1} : zero4f);

  auto mtype = 2;
  if (material.gltf_textures) mtype = 3;
  set_gluniform(state.program, "mat_type", mtype);
  set_gluniform(state.program, "mat_ke", material.emission);
  set_gluniform(state.program, "mat_kd", material.diffuse);
  set_gluniform(state.program, "mat_ks", vec3f{material.metallic});
  set_gluniform(state.program, "mat_rs", material.roughness);
  set_gluniform(state.program, "mat_op", material.opacity);
  set_gluniform(state.program, "mat_double_sided", (int)options.double_sided);
  set_gluniform_texture(state.program, "mat_ke_txt", "mat_ke_txt_on",
      material.emission_tex >= 0 ? state.textures.at(material.emission_tex)
                                 : opengl_texture{},
      0);
  set_gluniform_texture(state.program, "mat_kd_txt", "mat_kd_txt_on",
      material.diffuse_tex >= 0 ? state.textures.at(material.diffuse_tex)
                                : opengl_texture{},
      1);
  set_gluniform_texture(state.program, "mat_ks_txt", "mat_ks_txt_on",
      material.metallic_tex >= 0 ? state.textures.at(material.metallic_tex)
                                 : opengl_texture{},
      2);
  set_gluniform_texture(state.program, "mat_rs_txt", "mat_rs_txt_on",
      material.roughness_tex >= 0 ? state.textures.at(material.roughness_tex)
                                  : opengl_texture{},
      3);
  set_gluniform_texture(state.program, "mat_norm_txt", "mat_norm_txt_on",
      material.normal_tex >= 0 ? state.textures.at(material.normal_tex)
                               : opengl_texture{},
      5);

  set_gluniform(state.program, "elem_faceted", (int)shape.normals.empty());
  set_glvertexattrib(state.program, "vert_pos", vbos.positions_buffer, zero3f);
  set_glvertexattrib(state.program, "vert_norm", vbos.normals_buffer, zero3f);
  set_glvertexattrib(
      state.program, "vert_texcoord", vbos.texcoords_buffer, zero2f);
  set_glvertexattrib(
      state.program, "vert_color", vbos.colors_buffer, vec4f{1, 1, 1, 1});
  set_glvertexattrib(state.program, "vert_tangsp", vbos.tangentspaces_buffer,
      vec4f{0, 0, 1, 1});

  if (vbos.points_buffer) {
    set_gluniform(state.program, "elem_type", 1);
    draw_glpoints(vbos.points_buffer, vbos.points_buffer.num);
  }
  if (vbos.lines_buffer) {
    set_gluniform(state.program, "elem_type", 2);
    draw_gllines(vbos.lines_buffer, vbos.lines_buffer.num);
  }
  if (vbos.triangles_buffer) {
    set_gluniform(state.program, "elem_type", 3);
    draw_gltriangles(vbos.triangles_buffer, vbos.triangles_buffer.num);
  }
  if (vbos.quads_buffer) {
    set_gluniform(state.program, "elem_type", 3);
    draw_gltriangles(vbos.quads_buffer, vbos.quads_buffer.num);
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
  if (options.edges) throw std::runtime_error("edges are momentarily disabled");

  // for (int i = 0; i < 16; i++) { glDisableVertexAttribArray(i); }
}

// Display a scene
void draw_glscene(drawgl_state& state, const yocto_scene& scene,
    const vec4i& viewport, const app_selection& highlighted,
    const drawgl_params& options) {
  auto& camera      = scene.cameras.at(options.camera);
  auto  camera_view = mat4f(inverse(camera.frame));
  auto  camera_proj = perspective_mat(
      camera_fov(camera).x * (float)viewport.w / (float)viewport.z,
      (float)viewport.z / (float)viewport.w, options.near, options.far);

  bind_glprogram(state.program);
  set_gluniform(state.program, "cam_pos", camera.frame.o);
  set_gluniform(state.program, "cam_xform_inv", camera_view);
  set_gluniform(state.program, "cam_proj", camera_proj);
  set_gluniform(state.program, "eyelight", (int)options.eyelight);
  set_gluniform(state.program, "exposure", options.exposure);
  set_gluniform(state.program, "gamma", options.gamma);

  if (!options.eyelight) {
    auto lights_pos  = vector<vec3f>();
    auto lights_ke   = vector<vec3f>();
    auto lights_type = vector<int>();
    for (auto& instance : scene.instances) {
      if (instance.shape < 0) continue;
      auto& shape    = scene.shapes[instance.shape];
      auto& material = scene.materials[instance.material];
      if (material.emission == zero3f) continue;
      if (lights_pos.size() >= 16) break;
      auto bbox = compute_bounds(shape);
      auto pos  = (bbox.max + bbox.min) / 2;
      auto area = 0.0f;
      if (!shape.triangles.empty()) {
        for (auto t : shape.triangles)
          area += triangle_area(
              shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
      } else if (!shape.quads.empty()) {
        for (auto q : shape.quads)
          area += quad_area(shape.positions[q.x], shape.positions[q.y],
              shape.positions[q.z], shape.positions[q.w]);
      } else if (!shape.lines.empty()) {
        for (auto l : shape.lines)
          area += line_length(shape.positions[l.x], shape.positions[l.y]);
      } else {
        area += shape.positions.size();
      }
      auto ke = material.emission * area;
      lights_pos.push_back(transform_point(instance.frame, pos));
      lights_ke.push_back(ke);
      lights_type.push_back(0);
    }
    set_gluniform(state.program, "lamb", zero3f);
    set_gluniform(state.program, "lnum", (int)lights_pos.size());
    for (auto i = 0; i < lights_pos.size(); i++) {
      auto is = std::to_string(i);
      set_gluniform(state.program, ("lpos[" + is + "]").c_str(), lights_pos[i]);
      set_gluniform(state.program, ("lke[" + is + "]").c_str(), lights_ke[i]);
      set_gluniform(
          state.program, ("ltype[" + is + "]").c_str(), (int)lights_type[i]);
    }
  }

  if (options.wireframe) set_glwireframe(true);
  for (auto instance_id = 0; instance_id < scene.instances.size();
       instance_id++) {
    auto& instance = scene.instances[instance_id];
    // auto& shape     = scene.shapes[instance.shape];
    // auto& material  = scene.materials[shape.material];
    auto highlight = highlighted.type == typeid(yocto_instance) &&
                     highlighted.index == instance_id;
    draw_glinstance(state, scene, instance, highlight, options);
  }

  unbind_opengl_program();
  if (options.wireframe) set_glwireframe(false);
}

void init_drawgl_state(drawgl_state& state, const yocto_scene& scene) {
  // load textures and vbos
  init_glprogram(state.program, vertex, fragment);
  state.textures.resize(scene.textures.size());
  for (auto texture_id = 0; texture_id < scene.textures.size(); texture_id++) {
    auto& texture = scene.textures[texture_id];
    if (!texture.hdr.empty()) {
      init_gltexture(state.textures[texture_id], texture.hdr, true, true, true);
    } else if (!texture.ldr.empty()) {
      init_gltexture(state.textures[texture_id], texture.ldr, true, true, true);
    } else {
      throw std::runtime_error("bad texture");
    }
  }
  state.shapes.resize(scene.shapes.size());
  for (auto shape_id = 0; shape_id < scene.shapes.size(); shape_id++) {
    auto& shape = scene.shapes[shape_id];
    auto  vbos  = drawgl_shape();
    if (shape.quadspos.empty()) {
      if (!shape.positions.empty())
        init_glarraybuffer(vbos.positions_buffer, shape.positions, false);
      if (!shape.normals.empty())
        init_glarraybuffer(vbos.normals_buffer, shape.normals, false);
      if (!shape.texcoords.empty())
        init_glarraybuffer(vbos.texcoords_buffer, shape.texcoords, false);
      if (!shape.colors.empty())
        init_glarraybuffer(vbos.colors_buffer, shape.colors, false);
      if (!shape.tangents.empty())
        init_glarraybuffer(vbos.tangentspaces_buffer, shape.tangents, false);
      if (!shape.points.empty())
        init_glelementbuffer(vbos.points_buffer, shape.points, false);
      if (!shape.lines.empty())
        init_glelementbuffer(vbos.lines_buffer, shape.lines, false);
      if (!shape.triangles.empty())
        init_glelementbuffer(vbos.triangles_buffer, shape.triangles, false);
      if (!shape.quads.empty()) {
        auto triangles = quads_to_triangles(shape.quads);
        init_glelementbuffer(vbos.quads_buffer, triangles, false);
      }
    } else {
      auto quads     = vector<vec4i>{};
      auto positions = vector<vec3f>{};
      auto normals   = vector<vec3f>{};
      auto texcoords = vector<vec2f>{};
      split_facevarying(quads, positions, normals, texcoords, shape.quadspos,
          shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
          shape.texcoords);
      if (!positions.empty())
        init_glarraybuffer(vbos.positions_buffer, positions, false);
      if (!normals.empty())
        init_glarraybuffer(vbos.normals_buffer, normals, false);
      if (!texcoords.empty())
        init_glarraybuffer(vbos.texcoords_buffer, texcoords, false);
      if (!quads.empty()) {
        auto triangles = quads_to_triangles(quads);
        init_glelementbuffer(vbos.quads_buffer, triangles, false);
      }
    }
    state.shapes[shape_id] = vbos;
  }
}

void delete_drawgl_shape(drawgl_shape& glshape) {
  delete_glarraybuffer(glshape.positions_buffer);
  delete_glarraybuffer(glshape.normals_buffer);
  delete_glarraybuffer(glshape.texcoords_buffer);
  delete_glarraybuffer(glshape.colors_buffer);
  delete_glarraybuffer(glshape.tangentspaces_buffer);
  delete_glelementbuffer(glshape.points_buffer);
  delete_glelementbuffer(glshape.lines_buffer);
  delete_glelementbuffer(glshape.triangles_buffer);
  delete_glelementbuffer(glshape.quads_buffer);
}

// delete state
void delete_drawgl_state(drawgl_state& state) {
  if (!state.program) return;
  delete_glprogram(state.program);
  for (auto& texture : state.textures) delete_gltexture(texture);
  for (auto& shape : state.shapes) delete_drawgl_shape(shape);
  state.textures.clear();
  state.shapes.clear();
}

// draw with shading
void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app = *(app_state*)get_gluser_pointer(win);
  if (!begin_glwidgets_window(win, "yscnview")) return;
  if (!app.errors.empty() && error_message.empty()) {
    error_message = app.errors.front();
    app.errors.pop_front();
    open_glmodal(win, "error");
  }
  if (!draw_glmessage(win, "error", error_message)) {
    error_message = "";
  }
  if (draw_glfiledialog(
          win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
    add_new_scene(app, load_path);
  }
  if (draw_glfiledialog(win, "save", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.yaml;*.obj;*.pbrt")) {
    app.scenes[app.selected].outname = save_path;
    app.scenes[app.selected].task_queue.emplace_back(app_task_type::save_scene);
    save_path = "";
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    app.scenes[app.selected].imagename = save_path;
    app.scenes[app.selected].task_queue.emplace_back(app_task_type::save_image);
    save_path = "";
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save",
          app.selected >= 0 && app.scenes[app.selected].task_queue.empty())) {
    save_path = app.scenes[app.selected].outname;
    open_glmodal(win, "save");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save image", app.selected >= 0)) {
    save_path = app.scenes[app.selected].imagename;
    open_glmodal(win, "save image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", app.selected >= 0)) {
    app.scenes[app.selected].task_queue.emplace_back(
        app_task_type::close_scene);
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.scenes.empty()) return;
  draw_glcombobox(
      win, "scene", app.selected, (int)app.scenes.size(),
      [&app](int idx) { return app.scenes[idx].name.c_str(); }, false);
  auto& scn = app.scenes[app.selected];
  if (begin_glheader(win, "view")) {
    auto cam_names = vector<string>();
    for (auto& camera : scn.scene.cameras) cam_names.push_back(camera.uri);
    auto drawgl_prms = scn.drawgl_prms;
    if (scn.load_done) {
      if (draw_glcombobox(win, "camera", drawgl_prms.camera, cam_names)) {
      }
    }
    draw_glslider(win, "resolution", drawgl_prms.resolution, 0, 4096);
    draw_glcheckbox(win, "eyelight", drawgl_prms.eyelight);
    continue_glline(win);
    draw_glcheckbox(win, "wireframe", drawgl_prms.wireframe);
    continue_glline(win);
    draw_glcheckbox(win, "edges", drawgl_prms.edges);
    if (scn.time_range != zero2f) {
      draw_glslider(win, "time", scn.time, scn.time_range.x, scn.time_range.y);
      draw_gltextinput(win, "anim group", scn.anim_group);
      draw_glcheckbox(win, "animate", scn.animate);
    }
    draw_glslider(win, "exposure", drawgl_prms.exposure, -10, 10);
    draw_glslider(win, "gamma", drawgl_prms.gamma, 0.1f, 4);
    draw_glcheckbox(win, "double sided", drawgl_prms.double_sided);
    draw_glslider(win, "near", drawgl_prms.near, 0.01f, 1.0f);
    draw_glslider(win, "far", drawgl_prms.far, 1000.0f, 10000.0f);

    if (drawgl_prms != scn.drawgl_prms) {
      scn.task_queue.emplace_back(app_task_type::apply_edit,
          app_edit{typeid(drawgl_params), -1, drawgl_prms, false});
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    draw_gllabel(win, "scene", fs::path(scn.filename).filename());
    draw_gllabel(win, "filename", scn.filename);
    draw_gllabel(win, "outname", scn.outname);
    draw_gllabel(win, "imagename", scn.imagename);
    continue_glline(win);
    draw_glcheckbox(win, "fps", scn.navigation_fps);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : scn.scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      printf("%s\n", format_stats(scn.scene).c_str());
    }
    end_glheader(win);
  }
  if (scn.load_done && begin_glheader(win, "scene tree")) {
    draw_glscenetree(win, "", scn.scene, scn.selection, 200);
    end_glheader(win);
  }
  if (scn.load_done && begin_glheader(win, "scene object")) {
    auto edit = app_edit{};
    if (draw_glsceneinspector(win, "", scn.scene, scn.selection, edit, 200)) {
      scn.task_queue.emplace_back(app_task_type::apply_edit, edit);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);

  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  set_glviewport(get_glframebuffer_viewport(win));
  if (!app.scenes.empty() && app.selected >= 0 &&
      app.scenes[app.selected].load_done) {
    auto& scn = app.scenes[app.selected];
    draw_glscene(scn.state, scn.scene, get_glframebuffer_viewport(win),
        scn.selection, scn.drawgl_prms);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void apply_edit(const string& filename, yocto_scene& scene,
    drawgl_lights& lights, drawgl_params& drawgl_prms, float time,
    const string& anim_group, bool& reload_element, const app_edit& edit) {
  static auto last_time             = 0.0f;
  bool        updated_hierarchy     = false;
  auto& [type, index, data, reload] = edit;
  if (type == typeid(yocto_camera)) {
    scene.cameras[index] = any_cast<yocto_camera>(data);
  } else if (type == typeid(yocto_texture)) {
    scene.textures[index] = any_cast<yocto_texture>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_voltexture)) {
    scene.voltextures[index] = any_cast<yocto_voltexture>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_shape)) {
    scene.shapes[index] = any_cast<yocto_shape>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_subdiv)) {
    // TODO: this needs more fixing?
    scene.subdivs[index] = any_cast<yocto_subdiv>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_material)) {
    scene.materials[index] = any_cast<yocto_material>(data);
  } else if (type == typeid(yocto_instance)) {
    scene.instances[index] = any_cast<yocto_instance>(data);
  } else if (type == typeid(yocto_animation)) {
    scene.animations[index] = any_cast<yocto_animation>(data);
    updated_hierarchy       = true;
  } else if (type == typeid(yocto_scene_node)) {
    scene.nodes[index] = any_cast<yocto_scene_node>(data);
    updated_hierarchy  = true;
  } else if (type == typeid(drawgl_params)) {
    drawgl_prms = any_cast<drawgl_params>(data);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }
  if (updated_hierarchy || time != last_time) {
    update_transforms(scene, time, anim_group);
    last_time = time;
  }
}

// reload an element
void load_element(
    const string& filename, yocto_scene& scene, const app_edit& edit) {
  auto& [type, index, data, reload] = edit;

  if (type == typeid(yocto_texture)) {
    auto& texture = scene.textures[index];
    if (is_hdr_filename(texture.uri)) {
      load_image(fs::path(filename).parent_path() / texture.uri, texture.hdr);
    } else {
      load_imageb(fs::path(filename).parent_path() / texture.uri, texture.ldr);
    }
  } else if (type == typeid(yocto_voltexture)) {
    auto& texture = scene.voltextures[index];
    load_volume(fs::path(filename).parent_path() / texture.uri, texture.vol);
  } else if (type == typeid(yocto_shape)) {
    auto& shape = scene.shapes[index];
    load_shape(fs::path(filename).parent_path() / shape.uri, shape.points,
        shape.lines, shape.triangles, shape.quads, shape.quadspos,
        shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
        shape.texcoords, shape.colors, shape.radius, false);
  } else if (type == typeid(yocto_subdiv)) {
    // TODO: this needs more fixing?
    auto& subdiv = scene.subdivs[index];
    load_shape(fs::path(filename).parent_path() / subdiv.uri, subdiv.points,
        subdiv.lines, subdiv.triangles, subdiv.quads, subdiv.quadspos,
        subdiv.quadsnorm, subdiv.quadstexcoord, subdiv.positions,
        subdiv.normals, subdiv.texcoords, subdiv.colors, subdiv.radius,
        subdiv.facevarying);
    tesselate_subdiv(scene, scene.subdivs[index]);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }
}

// update
void update(const opengl_window& win, app_state& app) {
  // close if needed
  while (!app.scenes.empty()) {
    auto pos = -1;
    for (auto idx = 0; idx < app.scenes.size(); idx++) {
      for (auto& task : app.scenes[idx].task_queue) {
        if (task.type == app_task_type::close_scene) pos = idx;
      }
    }
    if (pos < 0) break;
    app.scenes.erase(app.scenes.begin() + pos);
    app.selected = app.scenes.empty() ? -1 : 0;
  }

  // apply synchronous edit
  for (auto& scn : app.scenes) {
    while (!scn.task_queue.empty()) {
      auto& task = scn.task_queue.front();
      if (task.type != app_task_type::apply_edit) break;
      log_glinfo(win, "start editing " + scn.filename);
      try {
        auto reload_element = false;
        apply_edit(scn.filename, scn.scene, scn.lights, scn.drawgl_prms,
            scn.time, scn.anim_group, reload_element, task.edit);
        log_glinfo(win, "done editing " + scn.filename);
        if (reload_element) {
          scn.load_done = false;
          scn.task_queue.emplace_back(app_task_type::load_element, task.edit);
        }
      } catch (std::exception& e) {
        log_glerror(win, e.what());
        app.errors.push_back("cannot edit " + scn.filename);
      }
      scn.task_queue.pop_front();
    }
  }

  // grab result of finished tasks
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (!task.result.valid() || task.result.wait_for(std::chrono::nanoseconds(
                                    10)) != std::future_status::ready)
      continue;
    switch (task.type) {
      case app_task_type::none: break;
      case app_task_type::close_scene: break;
      case app_task_type::load_scene: {
        try {
          task.result.get();
          scn.load_done = true;
          init_drawgl_state(scn.state, scn.scene);
          log_glinfo(win, "done loading " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [error]";
          app.errors.push_back("cannot load " + scn.filename);
        }
      } break;
      case app_task_type::load_element: {
        try {
          task.result.get();
          scn.load_done = true;
          log_glinfo(win, "done loading element from " + scn.filename);
          if (task.edit.type == typeid(yocto_texture)) {
            // not supported yet
            log_glerror(win, "texture refresh is not supported yet");
          } else if (task.edit.type == typeid(yocto_voltexture)) {
          } else if (task.edit.type == typeid(yocto_shape)) {
            // not supported yet
            log_glerror(win, "shape refresh is not supported yet");
          } else if (task.edit.type == typeid(yocto_subdiv)) {
            // not supported yet
            log_glerror(win, "shape refresh is not supported yet");
          } else {
            throw std::runtime_error("unsupported type");
          }
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [error]";
          app.errors.push_back("cannot load element from " + scn.filename);
        }
      } break;
      case app_task_type::save_image: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.imagename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.imagename);
        }
      } break;
      case app_task_type::save_scene: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.outname);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.outname);
        }
      } break;
      case app_task_type::apply_edit: break;
    }
    scn.task_queue.pop_front();
  }
  // schedule tasks not running
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (task.result.valid()) continue;
    task.stop = false;
    switch (task.type) {
      case app_task_type::none: break;
      case app_task_type::close_scene: break;
      case app_task_type::load_scene: {
        log_glinfo(win, "start loading " + scn.filename);
        scn.load_done = false;
        task.result   = std::async(std::launch::async, [&scn]() {
          load_scene(scn.filename, scn.scene, scn.load_prms);
          tesselate_subdivs(scn.scene);
          init_drawgl_lights(scn.lights, scn.scene);
          if (scn.lights.empty() && !scn.drawgl_prms.eyelight) {
            printf("no lights presents, switching to eyelight shader\n");
            scn.drawgl_prms.eyelight = true;
          }
          scn.time_range = compute_animation_range(scn.scene);
          scn.time       = scn.time_range.x;
        });
      } break;
      case app_task_type::load_element: {
        log_glinfo(win, "start loading element for " + scn.filename);
        scn.load_done = false;
        task.result   = std::async(std::launch::async, [&scn, &task]() {
          load_element(scn.filename, scn.scene, task.edit);
        });
      } break;
      case app_task_type::save_image: {
        log_glinfo(win, "start saving " + scn.imagename);
        task.result = std::async(std::launch::async,
            []() { throw std::runtime_error("not implemnted yet"); });
      } break;
      case app_task_type::save_scene: {
        log_glinfo(win, "start saving " + scn.outname);
        task.result = std::async(std::launch::async,
            [&scn]() { save_scene(scn.outname, scn.scene, scn.save_prms); });
      } break;
      case app_task_type::apply_edit: break;
    }
  }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  for (auto& path : paths) add_new_scene(app, path);
}

// run ui loop
void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnview", &app, draw);

  // init widget
  init_glwidgets(win);

  // loop
  auto mouse_pos = zero2f, last_pos = zero2f;
  auto last_time = std::chrono::high_resolution_clock::now();
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto alt_down       = get_glalt_key(win);
    auto shift_down     = get_glshift_key(win);
    auto widgets_active = get_glwidgets_active(win);

    // update trasforms
    if (app.selected >= 0) {
      auto& scn = app.scenes[app.selected];
      update_transforms(scn.scene, scn.time);
    }

    // handle mouse and keyboard for navigation
    if (app.selected >= 0 && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& scn        = app.scenes[app.selected];
      auto& old_camera = scn.scene.cameras.at(scn.drawgl_prms.camera);
      auto  camera     = scn.scene.cameras.at(scn.drawgl_prms.camera);
      auto  dolly      = 0.0f;
      auto  pan        = zero2f;
      auto  rotate     = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      if (camera.frame != old_camera.frame ||
          camera.focus != old_camera.focus) {
        scn.task_queue.emplace_back(app_task_type::apply_edit,
            app_edit{
                typeid(yocto_camera), scn.drawgl_prms.camera, camera, false});
      }
    }

    // animation
    if (app.selected >= 0 && app.scenes[app.selected].animate) {
      auto& scn     = app.scenes[app.selected];
      auto  now     = std::chrono::high_resolution_clock::now();
      auto  elapsed = now - last_time;
      auto  time    = (double)(elapsed.count()) / 1000000000.0;
      scn.time += min(1 / 60.0f, (float)time);
      if (scn.time < scn.time_range.x || scn.time > scn.time_range.y)
        scn.time = scn.time_range.x;
      update_transforms(scn.scene, scn.time);
      last_time = now;
    }

    // update
    update(win, app);

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // clear
  delete_glwindow(win);
}

int main(int argc, char* argv[]) {
  // initialize app
  app_state app{};
  auto      filenames  = vector<string>{};
  auto      noparallel = false;

  // parse command line
  auto parser = CLI::App{"views scenes inteactively"};
  parser.add_option("--camera", app.drawgl_prms.camera, "Camera index.");
  parser.add_option(
      "--resolution,-r", app.drawgl_prms.resolution, "Image resolution.");
  parser.add_flag("--eyelight,!--no-eyelight,-c", app.drawgl_prms.eyelight,
      "Eyelight rendering.");
  parser.add_flag("--noparallel", noparallel, "Disable parallel execution.");
  parser.add_option("scenes", filenames, "Scene filenames")->required(true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }

  // fix parallel code
  if (noparallel) {
    app.load_prms.noparallel = true;
    app.save_prms.noparallel = true;
  }

  // loading images
  for (auto filename : filenames) add_new_scene(app, filename);
  app.selected = app.scenes.empty() ? -1 : 0;

  // run ui
  run_ui(app);

  // done
  return 0;
}
