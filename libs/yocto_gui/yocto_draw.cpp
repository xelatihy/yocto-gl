//
// Utilities for real-time rendering of a scene.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include "yocto_draw.h"

#include <yocto/yocto_commonio.h>

#include <unordered_map>
#include <unordered_set>

#include "ext/glad/glad.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unordered_map;
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

namespace yocto {

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

static const char* glscene_vertex =
    R"(
#version 330

layout(location = 0) in vec3 positions;           // vertex position (in mesh coordinate frame)
layout(location = 1) in vec3 normals;             // vertex normal (in mesh coordinate frame)
layout(location = 2) in vec2 texcoords;           // vertex texcoords
layout(location = 3) in vec4 colors;              // vertex color
layout(location = 4) in vec4 tangents;            // vertex tangent space

uniform mat4 frame;             // shape transform
uniform mat4 frameit;           // shape transform
uniform float offset;           // shape normal offset

uniform mat4 view;              // inverse of the camera frame (as a matrix)
uniform mat4 projection;        // camera projection

out vec3 position;              // [to fragment shader] vertex position (in world coordinate)
out vec3 normal;                // [to fragment shader] vertex normal (in world coordinate)
out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
out vec4 color;                 // [to fragment shader] vertex color
out vec4 tangsp;                // [to fragment shader] vertex tangent space

// main function
void main() {
  // copy values
  position = positions;
  normal = normals;
  tangsp = tangents;
  texcoord = texcoords;
  color = colors;

  // normal offset
  if(offset != 0) {
    position += offset * normal;
  }

  // world projection
  position = (frame * vec4(position,1)).xyz;
  normal = (frameit * vec4(normal,0)).xyz;
  tangsp.xyz = (frame * vec4(tangsp.xyz,0)).xyz;

  // clip
  gl_Position = projection * view * vec4(position,1);
}
)";

static const char* glscene_fragment =
    R"(
#version 330

float pif = 3.14159265;

uniform bool eyelight;         // eyelight shading
uniform vec3 lamb;             // ambient light
uniform int  lnum;             // number of lights
uniform int  ltype[16];        // light type (0 -> point, 1 -> directional)
uniform vec3 lpos[16];         // light positions
uniform vec3 lke[16];          // light intensities

void evaluate_light(int lid, vec3 position, out vec3 cl, out vec3 wi) {
  cl = vec3(0,0,0);
  wi = vec3(0,0,0);
  if(ltype[lid] == 0) {
    // compute point light color at position
    cl = lke[lid] / pow(length(lpos[lid]-position),2);
    // compute light direction at position
    wi = normalize(lpos[lid]-position);
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
  } else if(etype == 3) {
      if(ndi<=0 || ndo <=0) return vec3(0);
      vec3 diff = ndi * kd / pif;
      if(ndh<=0) return diff;
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

uniform int etype;
uniform bool faceted;
uniform vec4 highlight;           // highlighted color

uniform int mtype;                // material type
uniform vec3 emission;            // material ke
uniform vec3 diffuse;             // material kd
uniform vec3 specular;            // material ks
uniform float roughness;          // material rs
uniform float opacity;            // material op

uniform bool emission_tex_on;     // material ke texture on
uniform sampler2D emission_tex;   // material ke texture
uniform bool diffuse_tex_on;      // material kd texture on
uniform sampler2D diffuse_tex;    // material kd texture
uniform bool specular_tex_on;     // material ks texture on
uniform sampler2D specular_tex;   // material ks texture
uniform bool roughness_tex_on;    // material rs texture on
uniform sampler2D roughness_tex;  // material rs texture
uniform bool opacity_tex_on;      // material op texture on
uniform sampler2D opacity_tex;    // material op texture

uniform bool mat_norm_tex_on;     // material normal texture on
uniform sampler2D mat_norm_tex;   // material normal texture

uniform bool double_sided;        // double sided rendering

uniform mat4 frame;              // shape transform
uniform mat4 frameit;            // shape transform

bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                    out vec3 kd, out vec3 ks, out float rs, out float op) {
  if(mtype == 0) {
    ke = emission;
    kd = vec3(0,0,0);
    ks = vec3(0,0,0);
    op = 1;
    return false;
  }

  ke = color.xyz * emission;
  kd = color.xyz * diffuse;
  ks = color.xyz * specular;
  rs = roughness;
  op = color.w * opacity;

  vec4 ke_tex = (emission_tex_on) ? texture(emission_tex,texcoord) : vec4(1,1,1,1);
  vec4 kd_tex = (diffuse_tex_on) ? texture(diffuse_tex,texcoord) : vec4(1,1,1,1);
  vec4 ks_tex = (specular_tex_on) ? texture(specular_tex,texcoord) : vec4(1,1,1,1);
  vec4 rs_tex = (roughness_tex_on) ? texture(roughness_tex,texcoord) : vec4(1,1,1,1);
  vec4 op_tex = (opacity_tex_on) ? texture(opacity_tex,texcoord) : vec4(1,1,1,1);

  // get material color from textures and adjust values
  ke *= ke_tex.xyz;
  vec3 kb = kd * kd_tex.xyz;
  float km = ks.x * ks_tex.z;
  kd = kb * (1 - km);
  ks = kb * km + vec3(0.04) * (1 - km);
  rs *= ks_tex.y;
  rs = rs*rs;
  op *= kd_tex.w;

  return true;
}

vec3 apply_normal_map(vec2 texcoord, vec3 normal, vec4 tangsp) {
    if(!mat_norm_tex_on) return normal;
  vec3 tangu = normalize((frame * vec4(normalize(tangsp.xyz),0)).xyz);
  vec3 tangv = normalize(cross(normal, tangu));
  if(tangsp.w < 0) tangv = -tangv;
  vec3 texture = 2 * pow(texture(mat_norm_tex,texcoord).xyz, vec3(1/2.2)) - 1;
  // texture.y = -texture.y;
  return normalize( tangu * texture.x + tangv * texture.y + normal * texture.z );
}

in vec3 position;              // [from vertex shader] position in world space
in vec3 normal;                // [from vertex shader] normal in world space (need normalization)
in vec2 texcoord;              // [from vertex shader] texcoord
in vec4 color;                 // [from vertex shader] color
in vec4 tangsp;                // [from vertex shader] tangent space

uniform vec3 eye;              // camera position
uniform mat4 view;             // inverse of the camera frame (as a matrix)
uniform mat4 projection;       // camera projection

uniform float exposure; 
uniform float gamma;

out vec4 frag_color;      

vec3 triangle_normal(vec3 position) {
  vec3 fdx = dFdx(position); 
  vec3 fdy = dFdy(position); 
  return normalize((frame * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

// main
void main() {
  // view vector
  vec3 wo = normalize(eye - position);

  // prepare normals
  vec3 n;
  if(faceted) {
    n = triangle_normal(position);
  } else {
    n = normalize(normal);
  }

  // apply normal map
  n = apply_normal_map(texcoord, n, tangsp);

  // use faceforward to ensure the normals points toward us
  if(double_sided) n = faceforward(n,-wo,n);

  // get material color from textures
  vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
  bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

  // exit if needed
  if(brdf_op < 0.005) discard;

  // check const color
  if(etype == 0) {
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
      c += pif * brdfcos((has_brdf) ? etype : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
    } else {
      // accumulate ambient
      c += lamb * brdf_kd;
      // foreach light
      for(int lid = 0; lid < lnum; lid ++) {
        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
        evaluate_light(lid, position, cl, wi);
        c += cl * brdfcos((has_brdf) ? etype : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
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

gui_scene::~gui_scene() {
  clear_scene(this);
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto texture : textures) delete texture;
  for (auto light : lights) delete light;
}

// Initialize an OpenGL scene
void init_scene(gui_scene* scene) {
  if (is_initialized(scene->program)) return;
  auto error = ""s, errorlog = ""s;
  init_program(
      scene->program, glscene_vertex, glscene_fragment, error, errorlog);
}
bool is_initialized(gui_scene* scene) { return is_initialized(scene->program); }

// Clear an OpenGL scene
void clear_scene(gui_scene* scene) {
  for (auto texture : scene->textures) clear_texture(texture);
  for (auto shape : scene->shapes) clear_shape(shape);
  clear_program(scene->program);
  clear_program(scene->environment_program);
  clear_cubemap(scene->environment_cubemap);
}

// add camera
gui_camera* add_camera(gui_scene* scene) {
  return scene->cameras.emplace_back(new gui_camera{});
}
void set_frame(gui_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(gui_camera* camera, float lens, float aspect, float film) {
  camera->lens   = lens;
  camera->aspect = aspect;
  camera->film   = film;
}
void set_nearfar(gui_camera* camera, float near, float far) {
  camera->near = near;
  camera->far  = far;
}

// add texture
ogl_texture* add_texture(gui_scene* scene) {
  return scene->textures.emplace_back(new ogl_texture{});
}

// add shape
ogl_shape* add_shape(gui_scene* scene) {
  auto shape = new ogl_shape{};
  set_shape(shape);
  scene->shapes.push_back(shape);
  return shape;
}

// add instance
gui_instance* add_instance(gui_scene* scene) {
  return scene->instances.emplace_back(new gui_instance{});
}
void set_frame(gui_instance* instance, const frame3f& frame) {
  instance->frame = frame;
}
void set_shape(gui_instance* instance, ogl_shape* shape) {
  instance->shape = shape;
}
void set_material(gui_instance* instance, gui_material* material) {
  instance->material = material;
}
void set_hidden(gui_instance* instance, bool hidden) {
  instance->hidden = hidden;
}
void set_highlighted(gui_instance* instance, bool highlighted) {
  instance->highlighted = highlighted;
}

// add material
gui_material* add_material(gui_scene* scene) {
  return scene->materials.emplace_back(new gui_material{});
}
void set_emission(
    gui_material* material, const vec3f& emission, ogl_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    gui_material* material, const vec3f& color, ogl_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    gui_material* material, float specular, ogl_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_roughness(
    gui_material* material, float roughness, ogl_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    gui_material* material, float opacity, ogl_texture* opacity_tex) {
  material->opacity = opacity;
}
void set_metallic(
    gui_material* material, float metallic, ogl_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_normalmap(gui_material* material, ogl_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// shortcuts
gui_camera* add_camera(gui_scene* scene, const frame3f& frame, float lens,
    float aspect, float film, float near, float far) {
  auto camera = add_camera(scene);
  set_frame(camera, frame);
  set_lens(camera, lens, aspect, film);
  set_nearfar(camera, near, far);
  return camera;
}
gui_material* add_material(gui_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    ogl_texture* emission_tex, ogl_texture* color_tex,
    ogl_texture* specular_tex, ogl_texture* metallic_tex,
    ogl_texture* roughness_tex, ogl_texture* normalmap_tex) {
  auto material = add_material(scene);
  set_emission(material, emission, emission_tex);
  set_color(material, color, color_tex);
  set_specular(material, specular, specular_tex);
  set_metallic(material, metallic, metallic_tex);
  set_roughness(material, roughness, roughness_tex);
  set_normalmap(material, normalmap_tex);
  return material;
}
ogl_shape* add_shape(gui_scene* scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, bool edges) {
  auto shape = add_shape(scene);
  set_points(shape, points);
  set_lines(shape, lines);
  set_triangles(shape, triangles);
  set_quads(shape, quads);
  set_positions(shape, positions);
  set_normals(shape, normals);
  set_texcoords(shape, texcoords);
  set_colors(shape, colors);
  if (edges && (!triangles.empty() || !quads.empty())) {
    set_edges(shape, triangles, quads);
  }
  return shape;
}
gui_instance* add_instance(gui_scene* scene, const frame3f& frame,
    ogl_shape* shape, gui_material* material, bool hidden, bool highlighted) {
  auto instance = add_instance(scene);
  set_frame(instance, frame);
  set_shape(instance, shape);
  set_material(instance, material);
  set_hidden(instance, hidden);
  set_highlighted(instance, highlighted);
  return instance;
}

// add light
gui_light* add_light(gui_scene* scene) {
  return scene->lights.emplace_back(new gui_light{});
}
void set_light(gui_light* light, const vec3f& position, const vec3f& emission,
    ogl_light_type type, bool camera) {
  light->position = position;
  light->emission = emission;
  light->type     = type;
  light->camera   = camera;
}
void clear_lights(gui_scene* scene) {
  for (auto light : scene->lights) delete light;
  scene->lights.clear();
}
bool has_max_lights(gui_scene* scene) { return scene->lights.size() >= 16; }
void add_default_lights(gui_scene* scene) {
  clear_lights(scene);
  set_light(add_light(scene), normalize(vec3f{1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{-1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{-1, -1, 1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true);
  set_light(add_light(scene), normalize(vec3f{0.1, 0.5, -1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true);
}

// Draw a shape
void draw_object(
    gui_scene* scene, gui_instance* instance, const gui_scene_params& params) {
  static auto empty_instances = vector<frame3f>{identity3x4f};

  if (instance->hidden) return;

  assert_ogl_error();
  auto shape_xform     = frame_to_mat(instance->frame);
  auto shape_inv_xform = transpose(
      frame_to_mat(inverse(instance->frame, params.non_rigid_frames)));
  set_uniform(scene->program, "frame", shape_xform);
  set_uniform(scene->program, "frameit", shape_inv_xform);
  set_uniform(scene->program, "offset", 0.0f);
  if (instance->highlighted) {
    set_uniform(scene->program, "highlight", vec4f{1, 1, 0, 1});
  } else {
    set_uniform(scene->program, "highlight", vec4f{0, 0, 0, 0});
  }
  assert_ogl_error();

  auto material = instance->material;
  auto mtype    = 2;
  set_uniform(scene->program, "mtype", mtype);
  set_uniform(scene->program, "emission", material->emission);
  set_uniform(scene->program, "diffuse", material->color);
  set_uniform(scene->program, "specular",
      vec3f{material->metallic, material->metallic, material->metallic});
  set_uniform(scene->program, "roughness", material->roughness);
  set_uniform(scene->program, "opacity", material->opacity);
  set_uniform(scene->program, "double_sided", (int)params.double_sided);
  set_uniform(scene->program, "emission_tex", "emission_tex_on",
      material->emission_tex, 0);
  set_uniform(
      scene->program, "diffuse_tex", "diffuse_tex_on", material->color_tex, 1);
  set_uniform(scene->program, "specular_tex", "specular_tex_on",
      material->metallic_tex, 2);
  set_uniform(scene->program, "roughness_tex", "roughness_tex_on",
      material->roughness_tex, 3);
  set_uniform(scene->program, "opacity_tex", "opacity_tex_on",
      material->opacity_tex, 4);
  set_uniform(scene->program, "mat_norm_tex", "mat_norm_tex_on",
      material->normal_tex, 5);
  assert_ogl_error();

  auto shape = instance->shape;
  // set_uniform(scene->program, "faceted", !is_initialized(shape->normals));
  // set_attribute(scene->program, "positions", shape->positions, vec3f{0, 0,
  // 0}); set_attribute(scene->program, "normals", shape->normals, vec3f{0, 0,
  // 1}); set_attribute(scene->program, "texcoords", shape->texcoords, vec2f{0,
  // 0}); set_attribute(scene->program, "colors", shape->colors, vec4f{1, 1, 1,
  // 1}); set_attribute(scene->program, "tangents", shape->tangents, vec4f{0, 0,
  // 1, 1});
  assert_ogl_error();

  if (is_initialized(shape->points)) {
    // glPointSize(shape->points_size);
    set_uniform(scene->program, "etype", 1);
    // draw_elements(shape->points);
  }
  if (is_initialized(shape->lines)) {
    set_uniform(scene->program, "etype", 2);
    // draw_elements(shape->lines);
  }
  if (is_initialized(shape->triangles)) {
    set_uniform(scene->program, "etype", 3);
    // draw_elements(shape->triangles);
  }
  if (is_initialized(shape->quads)) {
    set_uniform(scene->program, "etype", 3);
    // draw_elements(shape->quads);
  }
  draw_shape(shape);

  assert_ogl_error();

  if (is_initialized(shape->edges) && params.edges && !params.wireframe) {
    set_uniform(scene->program, "mtype", mtype);
    set_uniform(scene->program, "emission", vec3f{0, 0, 0});
    set_uniform(scene->program, "diffuse", vec3f{0, 0, 0});
    set_uniform(scene->program, "specular", vec3f{0, 0, 0});
    set_uniform(scene->program, "roughness", 0);
    set_uniform(scene->program, "etype", 3);
    // draw_elements(shape->edges);
    draw_shape(shape);
    assert_ogl_error();
  }
}

namespace ibl {
inline ogl_shape* cube_shape();
}

// Display a scene
void draw_scene(gui_scene* scene, gui_camera* camera, const vec4i& viewport,
    const gui_scene_params& params) {
  static auto camera_light0 = gui_light{normalize(vec3f{1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true};
  static auto camera_light1 = gui_light{normalize(vec3f{-1, 1, 1}),
      vec3f{pif / 2, pif / 2, pif / 2}, ogl_light_type::directional, true};
  static auto camera_light2 = gui_light{normalize(vec3f{-1, -1, 1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true};
  static auto camera_light3 = gui_light{normalize(vec3f{0.1, 0.5, -1}),
      vec3f{pif / 4, pif / 4, pif / 4}, ogl_light_type::directional, true};
  static auto camera_lights = vector<gui_light*>{
      &camera_light0, &camera_light1, &camera_light2, &camera_light3};
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * atan(camera->film / (camera_aspect * 2 * camera->lens)))
          : (2 * atan(camera->film / (2 * camera->lens)));
  auto camera_view = frame_to_mat(inverse(camera->frame));
  auto camera_proj = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  assert_ogl_error();
  clear_ogl_framebuffer(params.background);
  set_ogl_viewport(viewport);

  assert_ogl_error();
  bind_program(scene->program);
  assert_ogl_error();
  set_uniform(scene->program, "eye", camera->frame.o);
  set_uniform(scene->program, "view", camera_view);
  set_uniform(scene->program, "projection", camera_proj);
  set_uniform(scene->program, "eyelight",
      params.shading == gui_shading_type::eyelight ? 1 : 0);
  set_uniform(scene->program, "exposure", params.exposure);
  set_uniform(scene->program, "gamma", params.gamma);
  assert_ogl_error();

  if (params.shading == gui_shading_type::lights ||
      params.shading == gui_shading_type::camlights) {
    assert_ogl_error();
    auto& lights = params.shading == gui_shading_type::lights ? scene->lights
                                                              : camera_lights;
    set_uniform(scene->program, "lamb", vec3f{0, 0, 0});
    set_uniform(scene->program, "lnum", (int)lights.size());
    auto lid = 0;
    for (auto light : lights) {
      auto is = std::to_string(lid);
      if (light->camera) {
        auto position = light->type == ogl_light_type::directional
                            ? transform_direction(
                                  camera->frame, light->position)
                            : transform_point(camera->frame, light->position);
        set_uniform(scene->program, ("lpos[" + is + "]").c_str(), position);
      } else {
        set_uniform(
            scene->program, ("lpos[" + is + "]").c_str(), light->position);
      }
      set_uniform(scene->program, ("lke[" + is + "]").c_str(), light->emission);
      set_uniform(
          scene->program, ("ltype[" + is + "]").c_str(), (int)light->type);
      lid++;
    }
    assert_ogl_error();
  }

  if (params.wireframe) set_ogl_wireframe(true);
  for (auto instance : scene->instances) {
    draw_object(scene, instance, params);
  }

  bind_program(scene->environment_program);
  // set_scene_uniforms
  set_uniform(scene->environment_program, "eye", camera->frame.o);
  set_uniform(scene->environment_program, "view", camera_view);
  set_uniform(scene->environment_program, "projection", camera_proj);
  set_uniform(scene->environment_program, "exposure", params.exposure);
  set_uniform(scene->environment_program, "gamma", params.gamma);

  set_uniform(
      scene->environment_program, "environment", scene->environment_cubemap, 0);
  set_uniform(scene->environment_program, "roughness", 0.0f);

  auto cube = ibl::cube_shape();
  draw_shape(cube);
  // set_attribute(
  //     scene->environment_program, "positions", cube->positions, vec3f{0, 0,
  //     0});
  // draw_elements(cube->triangles);

  unbind_program();
  if (params.wireframe) set_ogl_wireframe(false);
}

// image based lighting
namespace ibl {

static ogl_program* load_program(
    const string& vertex_filename, const string& fragment_filename) {
  auto error           = ""s;
  auto vertex_source   = ""s;
  auto fragment_source = ""s;

  if (!load_text(vertex_filename, vertex_source, error)) {
    printf("error loading vertex shader (%s): \n%s\n", vertex_filename.c_str(),
        error.c_str());
    return nullptr;
  }
  if (!load_text(fragment_filename, fragment_source, error)) {
    printf("error loading fragment shader (%s): \n%s\n",
        fragment_filename.c_str(), error.c_str());
    return nullptr;
  }

  auto program   = new ogl_program();
  auto error_buf = ""s;
  if (!init_program(
          program, vertex_source, fragment_source, error, error_buf)) {
    printf("\nerror: %s\n", error.c_str());
    printf("    %s\n", error_buf.c_str());
    return nullptr;
  }
  return program;
}

inline ogl_shape* cube_shape() {
  // Do not call this function for the first time in a draw loop!
  // Untested behavior.
  static ogl_shape* cube = nullptr;
  if (!cube) {
    // clang-format off
    static const auto cube_positions = vector<vec3f>{
      {1, -1, -1}, {1, -1,  1}, {-1, -1,  1}, {-1, -1, -1},
      {1,  1, -1}, {1,  1,  1}, {-1,  1,  1}, {-1,  1, -1},
    };
    static const auto cube_triangles = vector<vec3i>{
      {1, 3, 0}, {7, 5, 4}, {4, 1, 0}, {5, 2, 1},
      {2, 7, 3}, {0, 7, 4}, {1, 2, 3}, {7, 6, 5},
      {4, 5, 1}, {5, 6, 2}, {2, 6, 7}, {0, 3, 7}
    };
    // clang-format on
    cube = new ogl_shape{};
    set_shape(cube);
    set_positions(cube, cube_positions);
    set_triangles(cube, cube_triangles);
  }
  return cube;
}

inline ogl_shape* brdf_plane_shape() {
  // Do not call this function for the first time in a draw loop!
  // Untested behavior.
  static ogl_shape* brdf_plane = nullptr;
  if (!brdf_plane) {
    // clang-format off
    static const auto brdf_plane_positions = vector<vec3f>{
      {-1, -1, 0}, {1, -1,  0}, {1, 1,  0}, {-1, 1, 0},
    };
    static const auto brdf_plane_triangles = vector<vec3i>{
      {0, 1, 3}, {3, 2, 1}
    };
    // clang-format on
    brdf_plane = new ogl_shape{};
    set_shape(brdf_plane);
    set_positions(brdf_plane, brdf_plane_positions);
    set_triangles(brdf_plane, brdf_plane_triangles);
  }
  return brdf_plane;
}

// Using 6 render passes, bake a cubemap given a sampler for the environment.
// The input sampler can be either a cubemap or a latlong texture.
template <typename Sampler>
inline void bake_cubemap(ogl_cubemap* cubemap, const Sampler* environment,
    ogl_program* program, int size, int num_mipmap_levels = 1) {
  // init cubemap with no data
  set_cubemap<float>(cubemap, size, 3, true, true, true);
  auto cube = cube_shape();

  auto framebuffer = ogl_framebuffer{};
  set_framebuffer(&framebuffer, {size, size});

  // clang-format off
  frame3f cameras[6] = {
    lookat_frame({0, 0, 0}, { 1, 0, 0}, {0, 1, 0}),
    lookat_frame({0, 0, 0}, {-1, 0, 0}, {0, 1, 0}),
    lookat_frame({0, 0, 0}, { 0,-1, 0}, {0, 0,-1}),
    lookat_frame({0, 0, 0}, { 0, 1, 0}, {0, 0, 1}),
    lookat_frame({0, 0, 0}, { 0, 0,-1}, {0, 1, 0}),
    lookat_frame({0, 0, 0}, { 0, 0, 1}, {0, 1, 0})
  };
  // clang-format on

  bind_framebuffer(&framebuffer);
  bind_program(program);
  for (int mipmap_level = 0; mipmap_level < num_mipmap_levels; mipmap_level++) {
    // resize render buffer and viewport
    set_framebuffer(&framebuffer, {size, size});
    set_ogl_viewport(vec2i{size, size});

    for (auto i = 0; i < 6; ++i) {
      // perspective_mat(fov, aspect, near, far)
      auto camera_proj = perspective_mat(radians(90), 1, 1, 100);
      auto camera_view = frame_to_mat(inverse(cameras[i]));

      set_framebuffer_texture(&framebuffer, cubemap, i, mipmap_level);
      clear_ogl_framebuffer({0, 0, 0, 0}, true);

      set_uniform(program, "view", camera_view);
      set_uniform(program, "projection", camera_proj);
      set_uniform(program, "eye", vec3f{0, 0, 0});
      set_uniform(program, "mipmap_level", mipmap_level);
      set_uniform(program, "environment", environment, 0);

      // set_attribute(program, "positions", cube->positions, vec3f{0, 0, 0});
      // draw_elements(cube->triangles);
      draw_shape(cube);
    }
    size /= 2;
  }
  unbind_program();
  unbind_framebuffer();
}

inline void bake_specular_brdf_texture(ogl_texture* texture) {
  auto size        = 512;
  auto framebuffer = ogl_framebuffer{};
  auto brdf_plane  = brdf_plane_shape();

  auto program = load_program(
      "apps/ibl/shaders/bake_brdf.vert", "apps/ibl/shaders/bake_brdf.frag");
  assert_ogl_error();

  // ************ create texture ***************
  texture->is_float  = true;
  texture->linear    = true;
  texture->nchannels = 3;
  texture->size      = {size, size};
  glGenTextures(1, &texture->texture_id);

  // pre-allocate enough memory for the LUT texture.
  glBindTexture(GL_TEXTURE_2D, texture->texture_id);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, size, size, 0, GL_RGB, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  // TODO(giacomo): mipmaps?

  assert_ogl_error();

  set_framebuffer(&framebuffer, {size, size});
  set_framebuffer_texture(&framebuffer, texture, 0);

  bind_framebuffer(&framebuffer);
  bind_program(program);

  set_ogl_viewport(vec2i{size, size});
  clear_ogl_framebuffer({0, 0, 0, 0}, true);

  // set_attribute(program, "positions", brdf_plane->positions, vec3f{0, 0, 0});
  // draw_elements(brdf_plane->triangles);
  draw_shape(brdf_plane);
  assert_ogl_error();

  unbind_program();
  unbind_framebuffer();
}

void init_ibl_data(gui_scene* scene, const ogl_texture* environment_texture) {
  scene->ibl_program = ibl::load_program(
      "apps/ibl/shaders/scene.vert", "apps/ibl/shaders/ibl.frag");
  scene->environment_program = ibl::load_program(
      "apps/ibl/shaders/environment.vert", "apps/ibl/shaders/environment.frag");

  // make cubemap from environment texture
  {
    auto size    = environment_texture->size.y;
    auto program = ibl::load_program("apps/ibl/shaders/bake_cubemap.vert",
        "apps/ibl/shaders/bake_environment.frag");
    bake_cubemap(
        scene->environment_cubemap, environment_texture, program, size);
  }

  // bake irradiance map
  {
    auto program = ibl::load_program("apps/ibl/shaders/bake_cubemap.vert",
        "apps/ibl/shaders/bake_irradiance.frag");
    bake_cubemap(
        scene->irradiance_map, scene->environment_cubemap, program, 64);
  }

  // bake specular map
  {
    auto program = ibl::load_program("apps/ibl/shaders/bake_cubemap.vert",
        "apps/ibl/shaders/bake_specular.frag");
    bake_cubemap(
        scene->prefiltered_map, scene->environment_cubemap, program, 256, 6);
  }

  bake_specular_brdf_texture(scene->brdf_lut);
}
}  // namespace ibl

}  // namespace yocto
