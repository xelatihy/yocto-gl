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
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <list>

#ifdef _WIN32
#undef near
#undef far
#endif

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

// Application state
struct app_state {
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

  // editing
  pair<string, int> selection = {"camera", 0};
};

// Application state
struct app_states {
  // data
  std::list<app_state>    states;
  int                     selected = -1;
  std::list<string>       errors;
  std::list<app_state>    loading;
  std::list<future<void>> load_workers;

  // get image
  app_state& get_selected() {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }
  const app_state& get_selected() const {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }

  // default options
  load_params   load_prms   = {};
  save_params   save_prms   = {};
  drawgl_params drawgl_prms = {};
};

void load_scene_async(app_states& apps, const string& filename) {
  auto& app       = apps.loading.emplace_back();
  app.filename    = filename;
  app.imagename   = replace_extension(filename, ".png");
  app.outname     = replace_extension(filename, ".edited.yaml");
  app.name        = get_filename(app.filename);
  app.load_prms   = app.load_prms;
  app.save_prms   = app.save_prms;
  app.drawgl_prms = app.drawgl_prms;
  apps.load_workers.push_back(run_async([&app]() {
    load_scene(app.filename, app.scene);
    tesselate_subdivs(app.scene);
    init_drawgl_lights(app.lights, app.scene);
    app.time_range = compute_animation_range(app.scene);
    app.time       = app.time_range.x;
  }));
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
    const vec4i& viewport, const pair<string, int>& highlighted,
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
    auto highlight = highlighted.first == "instance" &&
                     highlighted.second == instance_id;
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

bool draw_glwidgets_camera(const opengl_window& win, app_state& scene, int id) {
  auto& camera = scene.scene.cameras[id];
  auto  edited = 0;
  edited += (int)draw_gltextinput(win, "name", camera.name);
  edited += (int)draw_glslider(win, "frame.x", camera.frame.x, -1, 1);
  edited += (int)draw_glslider(win, "frame.y", camera.frame.y, -1, 1);
  edited += (int)draw_glslider(win, "frame.z", camera.frame.z, -1, 1);
  edited += (int)draw_glslider(win, "frame.o", camera.frame.o, -10, 10);
  edited += (int)draw_glcheckbox(win, "ortho", camera.orthographic);
  edited += (int)draw_glslider(win, "lens", camera.lens, 0.01f, 1);
  edited += (int)draw_glslider(win, "film", camera.film, 0.01f, 0.1f);
  edited += (int)draw_glslider(win, "focus", camera.focus, 0.01f, 1000);
  edited += (int)draw_glslider(win, "aperture", camera.aperture, 0, 5);
  auto from         = camera.frame.o,
       to           = camera.frame.o - camera.focus * camera.frame.z;
  auto from_changed = draw_glslider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_glslider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    camera.frame = lookat_frame(from, to, {0, 1, 0});
    camera.focus = length(from - to);
    edited += 1;
  }
  return edited;
}

/// Visit struct elements.
bool draw_glwidgets_texture(
    const opengl_window& win, app_state& scene, int id) {
  auto& texture      = scene.scene.textures[id];
  auto  old_filename = texture.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", texture.name);
  edited += draw_gltextinput(win, "filename", texture.filename);
  draw_gllabel(
      win, "hdr", "%d x %d", texture.hdr.size().x, texture.hdr.size().y);
  draw_gllabel(
      win, "ldr", "%d x %d", texture.ldr.size().x, texture.ldr.size().y);
  if (edited && old_filename != texture.filename) {
    try {
      if (is_hdr_filename(texture.filename)) {
        load_image(texture.filename, texture.hdr);
      } else {
        load_imageb(texture.filename, texture.ldr);
      }
    } catch (std::exception& e) {
      push_glmessage("cannot load " + texture.filename);
      log_glinfo(win, "cannot load " + texture.filename);
      log_glinfo(win, e.what());
    }
    // TODO: update lights
  }
  return edited;
}

bool draw_glwidgets_material(
    const opengl_window& win, app_state& scene, int id) {
  auto& material = scene.scene.materials[id];
  auto  edited   = 0;
  edited += draw_gltextinput(win, "name", material.name);
  edited += draw_glhdrcoloredit(win, "emission", material.emission);
  edited += draw_glcoloredit(win, "diffuse", material.diffuse);
  edited += draw_glcoloredit(win, "specular", material.specular);
  edited += draw_glslider(win, "metallic", material.metallic, 0, 1);
  edited += draw_glslider(win, "roughness", material.roughness, 0, 1);
  edited += draw_glcoloredit(win, "coat", material.coat);
  edited += draw_glcoloredit(win, "transmission", material.transmission);
  edited += draw_glcheckbox(win, "refract", material.refract);
  edited += draw_glcoloredit(win, "vol transmission", material.voltransmission);
  edited += draw_glcoloredit(win, "vol meanfreepath", material.volmeanfreepath);
  edited += draw_glcoloredit(win, "vol scatter", material.volscatter);
  edited += draw_glcoloredit(win, "vol emission", material.volemission);
  edited += draw_glslider(win, "vol scale", material.volscale, 0, 1);
  edited += draw_glslider(win, "vol anisotropy", material.volanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", material.opacity, 0, 1);

  edited += draw_glcombobox(
      win, "emission_tex", material.emission_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, scene.scene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, scene.scene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_state& scene, int id) {
  auto& shape        = scene.scene.shapes[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  draw_gllabel(win, "points", "%ld", shape.points.size());
  draw_gllabel(win, "lines", "%ld", shape.lines.size());
  draw_gllabel(win, "triangles", "%ld", shape.triangles.size());
  draw_gllabel(win, "quads", "%ld", shape.quads.size());
  draw_gllabel(win, "quads pos", "%ld", shape.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", shape.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", shape.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", shape.positions.size());
  draw_gllabel(win, "norm", "%ld", shape.normals.size());
  draw_gllabel(win, "texcoord", "%ld", shape.texcoords.size());
  draw_gllabel(win, "color", "%ld", shape.colors.size());
  draw_gllabel(win, "radius", "%ld", shape.radius.size());
  draw_gllabel(win, "tangsp", "%ld", shape.tangents.size());
  if (edited && old_filename != shape.filename) {
    try {
      load_shape(shape.filename, shape.points, shape.lines, shape.triangles,
          shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords, shape.colors,
          shape.radius, false);
    } catch (std::exception& e) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
      log_glinfo(win, e.what());
    }
    // TODO: update mesh state
    // TODO: update lights
  }
  return edited;
}

inline bool draw_glwidgets_subdiv(
    const opengl_window& win, app_state& scene, int id) {
  auto& shape        = scene.scene.subdivs[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  edited += draw_glslider(win, "subdivisions", shape.subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", shape.catmullclark);
  edited += draw_glcheckbox(win, "smooth", shape.smooth);
  edited += draw_glcheckbox(win, "facevarying", shape.facevarying);
  edited += draw_glcombobox(
      win, "shape", shape.shape, scene.scene.shapes, true);
  edited += draw_glcombobox(win, "displacement_tex", shape.displacement_tex,
      scene.scene.textures, true);
  edited += draw_glslider(win, "displacement", shape.displacement, 0, 1);
  draw_gllabel(win, "points", "%ld", shape.points.size());
  draw_gllabel(win, "lines", "%ld", shape.lines.size());
  draw_gllabel(win, "triangles", "%ld", shape.triangles.size());
  draw_gllabel(win, "quads", "%ld", shape.quads.size());
  draw_gllabel(win, "quads pos", "%ld", shape.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", shape.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", shape.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", shape.positions.size());
  draw_gllabel(win, "norm", "%ld", shape.normals.size());
  draw_gllabel(win, "texcoord", "%ld", shape.texcoords.size());
  draw_gllabel(win, "color", "%ld", shape.colors.size());
  draw_gllabel(win, "radius", "%ld", shape.radius.size());
  if (edited && old_filename != shape.filename) {
    try {
      load_shape(shape.filename, shape.points, shape.lines, shape.triangles,
          shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords, shape.colors,
          shape.radius, false);
    } catch (std::exception& e) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
      log_glinfo(win, e.what());
    }
    tesselate_subdiv(scene.scene, shape);
    // TODO: update mesh state
    // TODO: update lights
  }
  if (edited && old_filename == shape.filename) {
    tesselate_subdiv(scene.scene, shape);
    // TODO: update mesh state
    // TODO: update lights
  }
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window& win, app_state& scene, int id) {
  auto& instance     = scene.scene.instances[id];
  auto  old_instance = instance;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame[0]", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", instance.shape, scene.scene.shapes, true);
  edited += draw_glcombobox(
      win, "material", instance.material, scene.scene.materials, true);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, app_state& scene, int id) {
  auto& environment = scene.scene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame[0]", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      scene.scene.textures, true);
  if (edited) {
    // TODO: update lights
  }
  return edited;
}

// draw with shading
void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         apps     = *(app_states*)get_gluser_pointer(win);
  auto          scene_ok = !apps.states.empty() && apps.selected >= 0;
  if (!begin_glwidgets_window(win, "yscnview")) return;
  draw_glmessages(win);
  if (draw_glfiledialog_button(win, "load", true, "load", load_path, false,
          "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save", scene_ok, "save", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.yaml;*.obj;*.pbrt")) {
    auto& app   = apps.get_selected();
    app.outname = save_path;
    try {
      save_scene(app.outname, app.scene);
    } catch (std::exception& e) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", scene_ok)) {
    auto it = apps.states.begin();
    std::advance(it, apps.selected);
    apps.states.erase(it);
    apps.selected = apps.states.empty() ? -1 : 0;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (apps.states.empty()) return;
  draw_glcombobox(
      win, "scene", apps.selected, (int)apps.states.size(),
      [&apps](int idx) {
        auto it = apps.states.begin();
        std::advance(it, apps.selected);
        return it->name.c_str();
      },
      false);
  if (scene_ok && begin_glheader(win, "view")) {
    auto& app    = apps.get_selected();
    auto& params = app.drawgl_prms;
    draw_glcombobox(win, "camera", params.camera, app.scene.cameras);
    draw_glslider(win, "resolution", params.resolution, 0, 4096);
    draw_glcheckbox(win, "eyelight", params.eyelight);
    continue_glline(win);
    draw_glcheckbox(win, "wireframe", params.wireframe);
    continue_glline(win);
    draw_glcheckbox(win, "edges", params.edges);
    if (app.time_range != zero2f) {
      draw_glslider(win, "time", app.time, app.time_range.x, app.time_range.y);
      draw_gltextinput(win, "anim group", app.anim_group);
      draw_glcheckbox(win, "animate", app.animate);
    }
    draw_glslider(win, "exposure", params.exposure, -10, 10);
    draw_glslider(win, "gamma", params.gamma, 0.1f, 4);
    draw_glcheckbox(win, "double sided", params.double_sided);
    draw_glslider(win, "near", params.near, 0.01f, 1.0f);
    draw_glslider(win, "far", params.far, 1000.0f, 10000.0f);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    auto& app = apps.get_selected();
    draw_gllabel(win, "scene", get_filename(app.filename));
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    draw_gllabel(win, "imagename", app.imagename);
    continue_glline(win);
    draw_glcheckbox(win, "fps", app.navigation_fps);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : app.scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : format_stats(app.scene)) print_info(stat);
    }
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "edit")) {
    static auto labels = vector<string>{"camera", "shape", "environment",
        "instance", "materials", "textures", "subdivs"};
    auto&       app    = apps.get_selected();
    if (draw_glcombobox(win, "selection##1", app.selection.first, labels))
      app.selection.second = 0;
    if (app.selection.first == "camera") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.cameras);
      draw_glwidgets_camera(win, app, app.selection.second);
    } else if (app.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.textures);
      draw_glwidgets_texture(win, app, app.selection.second);
    } else if (app.selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.materials);
      draw_glwidgets_material(win, app, app.selection.second);
    } else if (app.selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.subdivs);
      draw_glwidgets_shape(win, app, app.selection.second);
    } else if (app.selection.first == "subdiv") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.subdivs);
      draw_glwidgets_subdiv(win, app, app.selection.second);
    } else if (app.selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.instances);
      draw_glwidgets_instance(win, app, app.selection.second);
    } else if (app.selection.first == "environment") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.environments);
      draw_glwidgets_environment(win, app, app.selection.second);
    }
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(const opengl_window& win) {
  auto& apps = *(app_states*)get_gluser_pointer(win);

  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  set_glviewport(get_glframebuffer_viewport(win));
  if (!apps.states.empty() && apps.selected >= 0) {
    auto& app = apps.get_selected();
    draw_glscene(app.state, app.scene, get_glframebuffer_viewport(win),
        app.selection, app.drawgl_prms);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

// update
void update(const opengl_window& win, app_states& apps) {
  while (!apps.load_workers.empty() && is_ready(apps.load_workers.front())) {
    try {
      apps.load_workers.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load scene " + apps.loading.front().filename);
      log_glinfo(win, "cannot load scene " + apps.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    apps.states.splice(apps.states.end(), apps.loading, apps.loading.begin());
    apps.load_workers.pop_front();
    init_drawgl_state(apps.states.back().state, apps.states.back().scene);
    if (apps.selected < 0) apps.selected = (int)apps.states.size() - 1;
  }
#if 0
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
#endif
}

// run ui loop
void run_ui(app_states& apps) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnview", &apps, draw);
  set_drop_glcallback(
      win, [](const opengl_window& win, const vector<string>& paths) {
        auto& app = *(app_states*)get_gluser_pointer(win);
        for (auto& path : paths) load_scene_async(app, path);
      });

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
    auto scene_ok       = !apps.states.empty() && apps.selected >= 0;

    // update trasforms
    if (scene_ok) {
      auto& app = apps.get_selected();
      update_transforms(app.scene, app.time);
    }

    // handle mouse and keyboard for navigation
    if (scene_ok && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& app    = apps.get_selected();
      auto& camera = app.scene.cameras.at(app.drawgl_prms.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
    }

    // animation
    if (scene_ok && apps.get_selected().animate) {
      auto& app     = apps.get_selected();
      auto  now     = std::chrono::high_resolution_clock::now();
      auto  elapsed = now - last_time;
      auto  time    = (double)(elapsed.count()) / 1000000000.0;
      app.time += min(1 / 60.0f, (float)time);
      if (app.time < app.time_range.x || app.time > app.time_range.y)
        app.time = app.time_range.x;
      update_transforms(app.scene, app.time);
      last_time = now;
    }

    // update
    update(win, apps);

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // clear
  delete_glwindow(win);
}

int main(int argc, const char* argv[]) {
  // initialize app
  app_states app{};
  auto       filenames  = vector<string>{};
  auto       noparallel = false;

  // parse command line
  auto cli = make_cli("yscnview", "views scenes inteactively");
  add_cli_option(cli, "--camera", app.drawgl_prms.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app.drawgl_prms.resolution, "Image resolution.");
  add_cli_option(cli, "--eyelight/--no-eyelight,-c", app.drawgl_prms.eyelight,
      "Eyelight rendering.");
  add_cli_option(
      cli, "--noparallel", noparallel, "Disable parallel execution.");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (noparallel) {
    app.load_prms.noparallel = true;
    app.save_prms.noparallel = true;
  }

  // loading images
  for (auto filename : filenames) load_scene_async(app, filename);

  // run ui
  run_ui(app);

  // done
  return 0;
}
