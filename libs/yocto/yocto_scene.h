//
// # Yocto/Scene: Scene representation
//
// Yocto/Scene define a simple scene representation, and related utilities,
// used to write other libraries in Yocto/GL.
// Yocto/Scene is implemented in `yocto_scene.h` and `yocto_scene.cpp`.
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

#ifndef _YOCTO_SCENE_H_
#define _YOCTO_SCENE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_image.h"
#include "yocto_math.h"

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct scene_bvh_node {
  bbox3f  bbox;
  int32_t start;
  int16_t num;
  int8_t  axis;
  bool    internal;
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// Application data is not stored explicitly.
struct scene_bvh {
  vector<scene_bvh_node> nodes      = {};
  vector<int>            primitives = {};
};

// Camera based on a simple lens model. The camera is placed using a frame.
// Camera projection is described in photographic terms. In particular,
// we specify film size (35mm by default), film aspect ration,
// the lens' focal length, the focus distance and the lens aperture.
// All values are in meters. Here are some common aspect ratios used in video
// and still photography.
// 3:2    on 35 mm:  0.036 x 0.024
// 16:9   on 35 mm:  0.036 x 0.02025 or 0.04267 x 0.024
// 2.35:1 on 35 mm:  0.036 x 0.01532 or 0.05640 x 0.024
// 2.39:1 on 35 mm:  0.036 x 0.01506 or 0.05736 x 0.024
// 2.4:1  on 35 mm:  0.036 x 0.015   or 0.05760 x 0.024 (approx. 2.39 : 1)
// To compute good apertures, one can use the F-stop number from photography
// and set the aperture to focal length over f-stop.
struct sceneio_camera {
  string  name         = "";
  frame3f frame        = identity3x4f;
  bool    orthographic = false;
  float   lens         = 0.050;
  float   film         = 0.036;
  float   aspect       = 1.500;
  float   focus        = 10000;
  float   aperture     = 0;
};

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
struct sceneio_texture {
  string       name = "";
  image<vec4f> hdr  = {};
  image<vec4b> ldr  = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct sceneio_material {
  // material data
  string name = "";

  // material
  vec3f emission     = {0, 0, 0};
  vec3f color        = {0, 0, 0};
  float specular     = 0;
  float roughness    = 0;
  float metallic     = 0;
  float ior          = 1.5;
  vec3f spectint     = {1, 1, 1};
  float coat         = 0;
  float transmission = 0;
  float translucency = 0;
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  bool  thin         = true;

  // textures
  sceneio_texture* emission_tex     = nullptr;
  sceneio_texture* color_tex        = nullptr;
  sceneio_texture* specular_tex     = nullptr;
  sceneio_texture* metallic_tex     = nullptr;
  sceneio_texture* roughness_tex    = nullptr;
  sceneio_texture* transmission_tex = nullptr;
  sceneio_texture* translucency_tex = nullptr;
  sceneio_texture* spectint_tex     = nullptr;
  sceneio_texture* scattering_tex   = nullptr;
  sceneio_texture* coat_tex         = nullptr;
  sceneio_texture* opacity_tex      = nullptr;
  sceneio_texture* normal_tex       = nullptr;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct sceneio_shape {
  // shape data
  string name = "";

  // primitives
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // face-varying primitives
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};

  // subdivision data [experimental]
  int  subdivisions = 0;
  bool catmullclark = true;
  bool smooth       = true;

  // displacement data [experimental]
  float            displacement     = 0;
  sceneio_texture* displacement_tex = nullptr;

  // computed properties
  scene_bvh* bvh = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif

  // element cdf for sampling
  vector<float> elements_cdf = {};

  // cleanup
  ~sceneio_shape();
};

// Object.
struct sceneio_instance {
  // instance data
  string            name     = "";
  frame3f           frame    = identity3x4f;
  sceneio_shape*    shape    = nullptr;
  sceneio_material* material = nullptr;
};

// Environment map.
struct sceneio_environment {
  string           name         = "";
  frame3f          frame        = identity3x4f;
  vec3f            emission     = {0, 0, 0};
  sceneio_texture* emission_tex = nullptr;

  // computed properties
  vector<float> texels_cdf = {};
};

// Scene lights used during rendering. These are created automatically.
struct scene_light {
  sceneio_instance*    instance    = nullptr;
  sceneio_environment* environment = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct sceneio_scene {
  // scene elements
  vector<sceneio_camera*>      cameras      = {};
  vector<sceneio_instance*>    instances    = {};
  vector<sceneio_environment*> environments = {};
  vector<sceneio_shape*>       shapes       = {};
  vector<sceneio_texture*>     textures     = {};
  vector<sceneio_material*>    materials    = {};

  // additional information
  string name      = "";
  string copyright = "";

  // computed properties
  vector<scene_light*> lights = {};
  scene_bvh*           bvh    = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif

  // cleanup
  ~sceneio_scene();
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// add element to a scene
sceneio_camera*      add_camera(sceneio_scene* scene, const string& name = "");
sceneio_environment* add_environment(
    sceneio_scene* scene, const string& name = "");
sceneio_instance* add_instance(sceneio_scene* scene, const string& name = "");
sceneio_material* add_material(sceneio_scene* scene, const string& name = "");
sceneio_shape*    add_shape(sceneio_scene* scene, const string& name = "");
sceneio_texture*  add_texture(sceneio_scene* scene, const string& name = "");
sceneio_instance* add_complete_instance(
    sceneio_scene* scene, const string& name = "");

// set camera properties
void set_frame(sceneio_camera* camera, const frame3f& frame);
void set_lens(sceneio_camera* camera, float lens, float aspect, float film,
    bool ortho = false);
void set_focus(sceneio_camera* camera, float aperture, float focus);

// set instance properties
void set_frame(sceneio_instance* instance, const frame3f& frame);
void set_material(sceneio_instance* instance, sceneio_material* material);
void set_shape(sceneio_instance* instance, sceneio_shape* shape);

// set texture properties
void set_texture(sceneio_texture* texture, const image<vec4b>& img);
void set_texture(sceneio_texture* texture, const image<vec4f>& img);

// set material properties
void set_emission(sceneio_material* material, const vec3f& emission,
    sceneio_texture* emission_tex = nullptr);
void set_color(sceneio_material* material, const vec3f& color,
    sceneio_texture* color_tex = nullptr);
void set_specular(sceneio_material* material, float specular = 1,
    sceneio_texture* specular_tex = nullptr);
void set_ior(sceneio_material* material, float ior);
void set_metallic(sceneio_material* material, float metallic,
    sceneio_texture* metallic_tex = nullptr);
void set_transmission(sceneio_material* material, float transmission, bool thin,
    float trdepth, sceneio_texture* transmission_tex = nullptr);
void set_translucency(sceneio_material* material, float translucency, bool thin,
    float trdepth, sceneio_texture* translucency_tex = nullptr);
void set_roughness(sceneio_material* material, float roughness,
    sceneio_texture* roughness_tex = nullptr);
void set_opacity(sceneio_material* material, float opacity,
    sceneio_texture* opacity_tex = nullptr);
void set_thin(sceneio_material* material, bool thin);
void set_scattering(sceneio_material* material, const vec3f& scattering,
    float scanisotropy, sceneio_texture* scattering_tex = nullptr);
void set_normalmap(sceneio_material* material, sceneio_texture* normal_tex);

// set shape properties
void set_points(sceneio_shape* shape, const vector<int>& points);
void set_lines(sceneio_shape* shape, const vector<vec2i>& lines);
void set_triangles(sceneio_shape* shape, const vector<vec3i>& triangles);
void set_quads(sceneio_shape* shape, const vector<vec4i>& quads);
void set_fvquads(sceneio_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord);
void set_positions(sceneio_shape* shape, const vector<vec3f>& positions);
void set_normals(sceneio_shape* shape, const vector<vec3f>& normals);
void set_texcoords(sceneio_shape* shape, const vector<vec2f>& texcoords);
void set_colors(sceneio_shape* shape, const vector<vec4f>& colors);
void set_radius(sceneio_shape* shape, const vector<float>& radius);
void set_tangents(sceneio_shape* shape, const vector<vec4f>& tangents);
void set_subdivision(sceneio_shape* shape, int subdivisions, bool catmullclark,
    bool smooth = true);
void set_displacement(sceneio_shape* shape, float displacement,
    sceneio_texture* displacement_tex);

// set environment properties
void set_frame(sceneio_environment* environment, const frame3f& frame);
void set_emission(sceneio_environment* environment, const vec3f& emission,
    sceneio_texture* emission_tex = nullptr);

// add missing elements
void add_cameras(sceneio_scene* scene);
void add_radius(sceneio_scene* scene, float radius = 0.001f);
void add_materials(sceneio_scene* scene);
void add_sky(sceneio_scene* scene, float sun_angle = pif / 4);

// Trim all unused memory
void trim_memory(sceneio_scene* scene);

// Clone a scene
void clone_scene(sceneio_scene* dest, const sceneio_scene* scene);

// compute scene bounds
bbox3f compute_bounds(const sceneio_scene* scene);

// get named camera or default if name is empty
sceneio_camera* get_camera(const sceneio_scene* scene, const string& name = "");

}  // namespace yocto

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera.
ray3f eval_camera(
    const sceneio_camera* camera, const vec2f& image_uv, const vec2f& lens_uv);

// Evaluates a texture
vec2i texture_size(const sceneio_texture* texture);
vec4f lookup_texture(const sceneio_texture* texture, const vec2i& ij,
    bool ldr_as_linear = false);
vec4f eval_texture(const sceneio_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// Evaluate instance properties
vec3f eval_position(
    const sceneio_instance* instance, int element, const vec2f& uv);
vec3f eval_element_normal(const sceneio_instance* instance, int element);
vec3f eval_normal(
    const sceneio_instance* instance, int element, const vec2f& uv);
vec2f eval_texcoord(
    const sceneio_instance* instance, int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
    const sceneio_instance* instance, int element);
vec3f eval_normalmap(
    const sceneio_instance* instance, int element, const vec2f& uv);
vec3f eval_shading_normal(const sceneio_instance* instance, int element,
    const vec2f& uv, const vec3f& outgoing);
vec4f eval_color(
    const sceneio_instance* instance, int element, const vec2f& uv);

// Environment
vec3f eval_environment(
    const sceneio_environment* environment, const vec3f& direction);
vec3f eval_environment(const sceneio_scene* scene, const vec3f& direction);

// Material sample
struct scene_material_sample {
  vec3f emission     = {0, 0, 0};
  vec3f color        = {0, 0, 0};
  float specular     = 0;
  float roughness    = 0;
  float metallic     = 0;
  float ior          = 1.5;
  vec3f spectint     = {1, 1, 1};
  float coat         = 0;
  float transmission = 0;
  float translucency = 0;
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  bool  thin         = true;
  vec3f normalmap    = {0, 0, 1};
};

// Evaluates material and textures
scene_material_sample eval_material(
    const sceneio_material* material, const vec2f& texcoord);

// Material Bsdf parameters
struct scene_bsdf {
  // brdf lobes
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f metal        = {0, 0, 0};
  vec3f coat         = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  vec3f translucency = {0, 0, 0};
  vec3f refraction   = {0, 0, 0};
  float roughness    = 0;
  float ior          = 1;
  vec3f meta         = {0, 0, 0};
  vec3f metak        = {0, 0, 0};
  // weights
  float diffuse_pdf      = 0;
  float specular_pdf     = 0;
  float metal_pdf        = 0;
  float coat_pdf         = 0;
  float transmission_pdf = 0;
  float translucency_pdf = 0;
  float refraction_pdf   = 0;
};

// Eval material to obtain emission, brdf and opacity.
vec3f eval_emission(const sceneio_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
// Eval material to obatain emission, brdf and opacity.
scene_bsdf eval_bsdf(const sceneio_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
float      eval_opacity(const sceneio_instance* instance, int element,
         const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
// check if a brdf is a delta
bool is_delta(const scene_bsdf& bsdf);

// Material volume parameters
struct scene_vsdf {
  vec3f density    = {0, 0, 0};
  vec3f scatter    = {0, 0, 0};
  float anisotropy = 0;
};

// check if we have a volume
bool has_volume(const sceneio_instance* instance);
// evaluate volume
scene_vsdf eval_vsdf(
    const sceneio_instance* instance, int element, const vec2f& uv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-SCENE INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// Strategy used to build the bvh
enum struct scene_bvh_type {
  default_,
  highquality,
  middle,
  balanced,
#ifdef YOCTO_EMBREE
  embree_default,
  embree_highquality,
  embree_compact  // only for copy interface
#endif
};

// Params for scene bvh build
struct scene_bvh_params {
  scene_bvh_type bvh        = scene_bvh_type::default_;
  bool           noparallel = false;
};

// Progress callback called when loading.
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Build the bvh acceleration structure.
void init_bvh(sceneio_scene* scene, const scene_bvh_params& params,
    progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(sceneio_scene*       scene,
    const vector<sceneio_instance*>& updated_objects,
    const vector<sceneio_shape*>&    updated_shapes,
    const scene_bvh_params&          params);

// Results of intersect functions that include hit flag, the instance id,
// the shape element id, the shape element uv and intersection distance.
// Results values are set only if hit is true.
struct scene_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
scene_intersection intersect_scene_bvh(const sceneio_scene* scene,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
scene_intersection intersect_instance_bvh(const sceneio_instance* instance,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

// Make Cornell Box scene
void make_cornellbox(sceneio_scene* scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Return scene statistics as list of strings.
vector<string> scene_stats(const sceneio_scene* scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
    const sceneio_scene* scene, bool notextures = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_shapes(sceneio_scene* scene, progress_callback progress_cb = {});
void tesselate_shape(sceneio_shape* shape);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARDS COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

using scene_model [[deprecated]]       = sceneio_scene;
using scene_camera [[deprecated]]      = sceneio_camera;
using scene_texture [[deprecated]]     = sceneio_texture;
using scene_material [[deprecated]]    = sceneio_material;
using scene_shape [[deprecated]]       = sceneio_shape;
using scene_instance [[deprecated]]    = sceneio_instance;
using scene_environment [[deprecated]] = sceneio_environment;

}  // namespace yocto

#endif
