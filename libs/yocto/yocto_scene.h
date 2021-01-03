//
// # Yocto/Scene: Scene reepresentation
//
// Yocto/Scene defines a simple scene representation, and related utilities,
// mostly geared towards scene creation and serialization. Scene serialization
// is implemented in Yocto/SceneIO.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Handles to refer to scene elements
inline const int invalid_handle = -1;
using element_handle            = int;
using camera_handle             = int;
using texture_handle            = int;
using material_handle           = int;
using shape_handle              = int;
using instance_handle           = int;
using environment_handle        = int;

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
  image<vec4f> hdr = {};
  image<vec4b> ldr = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct scene_material {
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
  texture_handle emission_tex     = invalid_handle;
  texture_handle color_tex        = invalid_handle;
  texture_handle specular_tex     = invalid_handle;
  texture_handle metallic_tex     = invalid_handle;
  texture_handle roughness_tex    = invalid_handle;
  texture_handle transmission_tex = invalid_handle;
  texture_handle translucency_tex = invalid_handle;
  texture_handle spectint_tex     = invalid_handle;
  texture_handle scattering_tex   = invalid_handle;
  texture_handle coat_tex         = invalid_handle;
  texture_handle opacity_tex      = invalid_handle;
  texture_handle normal_tex       = invalid_handle;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct scene_shape {
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
  float          displacement     = 0;
  texture_handle displacement_tex = invalid_handle;

  // element cdf for sampling
  vector<float> elements_cdf = {};
  // shape is assigned at creation
  int shape_id = -1;
};

// Object.
struct scene_instance {
  // instance data
  frame3f         frame    = identity3x4f;
  shape_handle    shape    = invalid_handle;
  material_handle material = invalid_handle;
};

// Environment map.
struct scene_environment {
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {0, 0, 0};
  texture_handle emission_tex = invalid_handle;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene_scene {
  // additional information
  string name      = "";
  string copyright = "";

  // scene elements
  vector<scene_camera>      cameras      = {};
  vector<scene_instance>    instances    = {};
  vector<scene_environment> environments = {};
  vector<scene_shape>       shapes       = {};
  vector<scene_texture>     textures     = {};
  vector<scene_material>    materials    = {};

  // names (this will be cleanup significantly later)
  vector<string> camera_names      = {};
  vector<string> texture_names     = {};
  vector<string> material_names    = {};
  vector<string> shape_names       = {};
  vector<string> instance_names    = {};
  vector<string> environment_names = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// add element to a scene
camera_handle      add_camera(scene_scene& scene, const string& name = "");
environment_handle add_environment(scene_scene& scene, const string& name = "");
instance_handle    add_instance(scene_scene& scene, const string& name = "");
material_handle    add_material(scene_scene& scene, const string& name = "");
material_handle    add_shape(scene_scene& scene, const string& name = "");
texture_handle     add_texture(scene_scene& scene, const string& name = "");
instance_handle add_complete_instance(scene_scene& scene, const string& name);

// get element from a scene
scene_camera&      get_camera(scene_scene& scene, camera_handle handle);
scene_environment& get_environment(
    scene_scene& scene, environment_handle handle);
scene_instance& get_instance(scene_scene& scene, instance_handle handle);
scene_material& get_material(scene_scene& scene, material_handle handle);
scene_shape&    get_shape(scene_scene& scene, shape_handle handle);
scene_texture&  get_texture(scene_scene& scene, texture_handle handle);
scene_instance& get_complete_instance(
    scene_scene& scene, instance_handle handle);

// get element from a scene
const scene_camera& get_camera(const scene_scene& scene, camera_handle handle);
const scene_environment& get_environment(
    const scene_scene& scene, environment_handle handle);
const scene_instance& get_instance(
    const scene_scene& scene, instance_handle handle);
const scene_material& get_material(
    const scene_scene& scene, material_handle handle);
const scene_shape&   get_shape(const scene_scene& scene, shape_handle handle);
const scene_texture& get_texture(
    const scene_scene& scene, texture_handle handle);
const scene_instance& get_complete_instance(
    const scene_scene& scene, instance_handle handle);

// add missing elements
void add_cameras(scene_scene& scene);
void add_radius(scene_scene& scene, float radius = 0.001f);
void add_materials(scene_scene& scene);
void add_sky(scene_scene& scene, float sun_angle = pif / 4);

// Trim all unused memory
void trim_memory(scene_scene& scene);

// compute scene bounds
bbox3f compute_bounds(const scene_scene& scene);

// get named camera or default if name is empty
scene_camera&       get_camera(scene_scene& scene, const string& name);
const scene_camera& get_camera(const scene_scene& scene, const string& name);
camera_handle get_camera_handle(const scene_scene& scene, const string& name);

// get name
string get_camera_name(const scene_scene& scene, int idx);
string get_environment_name(const scene_scene& scene, int idx);
string get_shape_name(const scene_scene& scene, int idx);
string get_texture_name(const scene_scene& scene, int idx);
string get_instance_name(const scene_scene& scene, int idx);
string get_material_name(const scene_scene& scene, int idx);

string get_camera_name(const scene_scene& scene, const scene_camera& camera);
string get_environment_name(
    const scene_scene& scene, const scene_environment& environment);
string get_shape_name(const scene_scene& scene, const scene_shape& shape);
string get_texture_name(const scene_scene& scene, const scene_texture& texture);
string get_instance_name(
    const scene_scene& scene, const scene_instance& instance);
string get_material_name(
    const scene_scene& scene, const scene_material& material);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

// Progress callback called when loading.
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Apply subdivision and displacement rules.
void tesselate_shapes(
    scene_scene& scene, const progress_callback& progress_cb = {});
void tesselate_shape(scene_scene& scene, scene_shape& shape);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera.
ray3f eval_camera(
    const scene_camera& camera, const vec2f& image_uv, const vec2f& lens_uv);

// Evaluates a texture
vec2i texture_size(const scene_texture& texture);
vec4f lookup_texture(
    const scene_texture& texture, const vec2i& ij, bool ldr_as_linear = false);
vec4f eval_texture(const scene_texture& texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);
vec4f eval_texture(const scene_scene& scene, texture_handle texture,
    const vec2f& uv, bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// Evaluate instance properties
vec3f eval_position(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec3f eval_element_normal(
    const scene_scene& scene, const scene_instance& instance, int element);
vec3f eval_normal(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec2f eval_texcoord(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
    const scene_scene& scene, const scene_instance& instance, int element);
vec3f eval_normalmap(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec3f eval_shading_normal(const scene_scene& scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing);
vec4f eval_color(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv);

// Environment
vec3f eval_environment(const scene_scene& scene,
    const scene_environment& environment, const vec3f& direction);
vec3f eval_environment(const scene_scene& scene, const vec3f& direction);

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
scene_material_sample eval_material(const scene_scene& scene,
    const scene_material& material, const vec2f& texcoord);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

// Make Cornell Box scene
void make_cornellbox(scene_scene& scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Return scene statistics as list of strings.
vector<string> scene_stats(const scene_scene& scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
    const scene_scene& scene, bool notextures = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARDS COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

using scene_model [[deprecated]] = scene_scene;

using sceneio_scene       = scene_scene;
using sceneio_camera      = scene_camera;
using sceneio_texture     = scene_texture;
using sceneio_material    = scene_material;
using sceneio_shape       = scene_shape;
using sceneio_instance    = scene_instance;
using sceneio_environment = scene_environment;

}  // namespace yocto

#endif
