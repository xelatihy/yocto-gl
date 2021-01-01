//
// # Yocto/SceneIO: Scene serialization
//
// Yocto/Scene defines a simple scene representation, and related utilities,
// mostly geared towards scene creation and serialization.
// Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
// and a custom Json format.
// Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`,
// and depends on `cgltf.h`.
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

#ifndef _YOCTO_SCENEIO_H_
#define _YOCTO_SCENEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <functional>
#include <memory>
#include <string>
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
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

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

  // element cdf for sampling
  vector<float> elements_cdf = {};
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

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct sceneio_scene {
  // additional information
  string name      = "";
  string copyright = "";

  // scene elements
  vector<sceneio_camera*>      cameras      = {};
  vector<sceneio_instance*>    instances    = {};
  vector<sceneio_environment*> environments = {};
  vector<sceneio_shape*>       shapes       = {};
  vector<sceneio_texture*>     textures     = {};
  vector<sceneio_material*>    materials    = {};

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
    sceneio_scene* scene, const string& name);

// add missing elements
void add_cameras(sceneio_scene* scene);
void add_radius(sceneio_scene* scene, float radius = 0.001f);
void add_materials(sceneio_scene* scene);
void add_sky(sceneio_scene* scene, float sun_angle = pif / 4);

// Trim all unused memory
void trim_memory(sceneio_scene* scene);

// compute scene bounds
bbox3f compute_bounds(const sceneio_scene* scene);

// get named camera or default if name is empty
sceneio_camera* get_camera(const sceneio_scene* scene, const string& name = "");

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Progress callback called when loading.
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Load/save a scene in the supported formats. Throws on error.
// Calls the progress callback, if defined, as we process more data.
bool load_scene(const string& filename, sceneio_scene* scene, string& error,
    const progress_callback& progress_cb = {}, bool noparallel = false);
bool save_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb = {},
    bool noparallel = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_shapes(
    sceneio_scene* scene, const progress_callback& progress_cb = {});
void tesselate_shape(sceneio_shape* shape);

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
