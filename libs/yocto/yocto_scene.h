//
// # Yocto/Scene: Tiny library for a simple scene representation
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

#include <functional>
#include <memory>

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
  bbox3f bbox;
  int    start;
  short  num;
  bool   internal;
  byte   axis;
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
struct scene_camera {
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
struct scene_texture {
  string       name    = "";
  image<vec3f> colorf  = {};
  image<vec3b> colorb  = {};
  image<float> scalarf = {};
  image<byte>  scalarb = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct scene_material {
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
  scene_texture* emission_tex     = nullptr;
  scene_texture* color_tex        = nullptr;
  scene_texture* specular_tex     = nullptr;
  scene_texture* metallic_tex     = nullptr;
  scene_texture* roughness_tex    = nullptr;
  scene_texture* transmission_tex = nullptr;
  scene_texture* translucency_tex = nullptr;
  scene_texture* spectint_tex     = nullptr;
  scene_texture* scattering_tex   = nullptr;
  scene_texture* coat_tex         = nullptr;
  scene_texture* opacity_tex      = nullptr;
  scene_texture* normal_tex       = nullptr;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct scene_shape {
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
  vector<vec3f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};

  // subdivision data [experimental]
  int  subdivisions = 0;
  bool catmullclark = true;
  bool smooth       = true;

  // displacement data [experimental]
  float          displacement     = 0;
  scene_texture* displacement_tex = nullptr;

  // computed properties
  scene_bvh* bvh = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif

  // element cdf for sampling
  vector<float> elements_cdf = {};

  // cleanup
  ~scene_shape();
};

// Object.
struct scene_object {
  // object data
  string          name     = "";
  frame3f         frame    = identity3x4f;
  scene_shape*    shape    = nullptr;
  scene_material* material = nullptr;
};

// Environment map.
struct scene_environment {
  string         name         = "";
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {0, 0, 0};
  scene_texture* emission_tex = nullptr;

  // computed properties
  vector<float> texels_cdf = {};
};

// Scene lights used during rendering. These are created automatically.
struct scene_light {
  scene_object*      object      = nullptr;
  scene_environment* environment = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene_model {
  // scene elements
  vector<scene_camera*>      cameras      = {};
  vector<scene_object*>      objects      = {};
  vector<scene_environment*> environments = {};
  vector<scene_shape*>       shapes       = {};
  vector<scene_texture*>     textures     = {};
  vector<scene_material*>    materials    = {};

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
  ~scene_model();
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// add element to a scene
scene_camera*      add_camera(scene_model* scene, const string& name = "");
scene_environment* add_environment(scene_model* scene, const string& name = "");
scene_object*      add_object(scene_model* scene, const string& name = "");
scene_material*    add_material(scene_model* scene, const string& name = "");
scene_shape*       add_shape(scene_model* scene, const string& name = "");
scene_texture*     add_texture(scene_model* scene, const string& name = "");
scene_object* add_complete_object(scene_model* scene, const string& name = "");

// camera properties
void set_frame(scene_camera* camera, const frame3f& frame);
void set_lens(scene_camera* camera, float lens, float aspect, float film,
    bool ortho = false);
void set_focus(scene_camera* camera, float aperture, float focus);

// object properties
void set_frame(scene_object* object, const frame3f& frame);
void set_material(scene_object* object, scene_material* material);
void set_shape(scene_object* object, scene_shape* shape);

// texture properties
void set_texture(scene_texture* texture, const image<vec3b>& img);
void set_texture(scene_texture* texture, const image<vec3f>& img);
void set_texture(scene_texture* texture, const image<byte>& img);
void set_texture(scene_texture* texture, const image<float>& img);

// material properties
void set_emission(scene_material* material, const vec3f& emission,
    scene_texture* emission_tex = nullptr);
void set_color(scene_material* material, const vec3f& color,
    scene_texture* color_tex = nullptr);
void set_specular(scene_material* material, float specular = 1,
    scene_texture* specular_tex = nullptr);
void set_ior(scene_material* material, float ior);
void set_metallic(scene_material* material, float metallic,
    scene_texture* metallic_tex = nullptr);
void set_transmission(scene_material* material, float transmission, bool thin,
    float trdepth, scene_texture* transmission_tex = nullptr);
void set_translucency(scene_material* material, float translucency, bool thin,
    float trdepth, scene_texture* translucency_tex = nullptr);
void set_roughness(scene_material* material, float roughness,
    scene_texture* roughness_tex = nullptr);
void set_opacity(scene_material* material, float opacity,
    scene_texture* opacity_tex = nullptr);
void set_thin(scene_material* material, bool thin);
void set_scattering(scene_material* material, const vec3f& scattering,
    float scanisotropy, scene_texture* scattering_tex = nullptr);
void set_normalmap(scene_material* material, scene_texture* normal_tex);

// shape properties
void set_points(scene_shape* shape, const vector<int>& points);
void set_lines(scene_shape* shape, const vector<vec2i>& lines);
void set_triangles(scene_shape* shape, const vector<vec3i>& triangles);
void set_quads(scene_shape* shape, const vector<vec4i>& quads);
void set_positions(scene_shape* shape, const vector<vec3f>& positions);
void set_normals(scene_shape* shape, const vector<vec3f>& normals);
void set_texcoords(scene_shape* shape, const vector<vec2f>& texcoords);
void set_colors(scene_shape* shape, const vector<vec3f>& colors);
void set_radius(scene_shape* shape, const vector<float>& radius);
void set_tangents(scene_shape* shape, const vector<vec4f>& tangents);

// environment properties
void set_frame(scene_environment* environment, const frame3f& frame);
void set_emission(scene_environment* environment, const vec3f& emission,
    scene_texture* emission_tex = nullptr);

// add missing elements
void add_cameras(scene_model* scene);
void add_radius(scene_model* scene, float radius = 0.001f);
void add_materials(scene_model* scene);
void add_sky(scene_model* scene, float sun_angle = pif / 4);

// Trim all unused memory
void trim_memory(scene_model* scene);

// Clone a scene
void clone_scene(scene_model* dest, const scene_model* scene);

// compute scene bounds
bbox3f compute_bounds(const scene_model* scene);

// get named camera or default if name is empty
scene_camera* get_camera(const scene_model* scene, const string& name = "");

}  // namespace yocto

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera.
ray3f eval_camera(
    const scene_camera* camera, const vec2f& image_uv, const vec2f& lens_uv);

// Evaluates a texture
vec2i texture_size(const scene_texture* texture);
vec3f lookup_texture(
    const scene_texture* texture, const vec2i& ij, bool ldr_as_linear = false);
vec3f eval_texture(const scene_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// Evaluate object properties
vec3f eval_position(const scene_object* object, int element, const vec2f& uv);
vec3f eval_element_normal(const scene_object* object, int element);
vec3f eval_normal(const scene_object* object, int element, const vec2f& uv);
vec2f eval_texcoord(const scene_object* object, int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
    const scene_object* object, int element);
vec3f eval_normalmap(const scene_object* object, int element, const vec2f& uv);
vec3f eval_shading_normal(const scene_object* object, int element,
    const vec2f& uv, const vec3f& outgoing);
vec3f eval_color(const scene_object* object, int element, const vec2f& uv);

// Environment
vec3f eval_environment(
    const scene_environment* environment, const vec3f& direction);
vec3f eval_environment(const scene_model* scene, const vec3f& direction);

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
    const scene_material* material, const vec2f& texcoord);

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
vec3f eval_emission(const scene_object* object, int element, const vec2f& uv,
    const vec3f& normal, const vec3f& outgoing);
// Eval material to obatain emission, brdf and opacity.
scene_bsdf eval_bsdf(const scene_object* object, int element, const vec2f& uv,
    const vec3f& normal, const vec3f& outgoing);
float eval_opacity(const scene_object* object, int element, const vec2f& uv,
    const vec3f& normal, const vec3f& outgoing);
// check if a brdf is a delta
bool is_delta(const scene_bsdf& bsdf);

// Material volume parameters
struct scene_vsdf {
  vec3f density    = {0, 0, 0};
  vec3f scatter    = {0, 0, 0};
  float anisotropy = 0;
};

// check if we have a volume
bool has_volume(const scene_object* object);
// evaluate volume
scene_vsdf eval_vsdf(const scene_object* object, int element, const vec2f& uv);

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
void init_bvh(scene_model* scene, const scene_bvh_params& params,
    progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(scene_model*     scene,
    const vector<scene_object*>& updated_objects,
    const vector<scene_shape*>& updated_shapes, const scene_bvh_params& params);

// Results of intersect functions that include hit flag, the instance id,
// the shape element id, the shape element uv and intersection distance.
// Results values are set only if hit is true.
struct scene_intersection {
  int   object   = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
scene_intersection intersect_scene_bvh(const scene_model* scene,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
scene_intersection intersect_instance_bvh(const scene_model* object,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
scene_intersection intersect_instance_bvh(const scene_object* object,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

// Make Cornell Box scene
void make_cornellbox(scene_model* scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Return scene statistics as list of strings.
vector<string> scene_stats(const scene_model* scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
    const scene_model* scene, bool notextures = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_shapes(scene_model* scene, progress_callback progress_cb = {});
void tesselate_shape(scene_shape* shape);

}  // namespace yocto

#endif
