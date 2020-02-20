//
// # Yocto/SceneIO: Tiny library for Yocto/Scene input and output
//
// Yocto/SceneIO provides loading and saving functionality for scenes
// in Yocto/GL. We support a simple to use JSON format, PLY, OBJ and glTF.
// The JSON serialization is a straight copy of the in-memory scene data.
// To speed up testing, we also support a binary format that is a dump of
// the current scene. This format should not be use for archival though.
//
//
// ## Scene Loading and Saving
//
// 1. load a scene with `load_scene()` and save it with `save_scene()`
//
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

#ifndef _YOCTO_SCENEIO_H_
#define _YOCTO_SCENEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <functional>
#include <memory>

#include "yocto_image.h"
#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Using directives.
using std::string;
using std::vector;
using namespace std::string_literals;
using std::function;
using yocto::image::image;
using namespace yocto::math;

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
struct camera {
  string  name         = "";
  frame3f frame        = identity3x4f;
  bool    orthographic = false;
  float   lens         = 0.050;
  float   film         = 0.036;
  float   aspect       = 1.500;
  float   focus        = flt_max;
  float   aperture     = 0;
};

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
struct texture {
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
struct material {
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
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  float displacement = 0;
  bool  thin         = true;

  // textures
  texture* emission_tex     = nullptr;
  texture* color_tex        = nullptr;
  texture* specular_tex     = nullptr;
  texture* metallic_tex     = nullptr;
  texture* roughness_tex    = nullptr;
  texture* transmission_tex = nullptr;
  texture* spectint_tex     = nullptr;
  texture* scattering_tex   = nullptr;
  texture* coat_tex         = nullptr;
  texture* opacity_tex      = nullptr;
  texture* normal_tex       = nullptr;
  texture* displacement_tex = nullptr;

  // [experimental] properties to drive subdiv and displacement
  int  subdivisions = 2;
  bool smooth       = true;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct shape {
  // shape data
  string name = "";

  // primitives
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec3f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};
};

// Subdiv data represented as indexed meshes of elements.
// May contain points, lines, triangles, quads or
// face-varying quads.
struct subdiv {
  // shape data
  string name = "";

  // face-varying primitives
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
};

// Instance data.
struct instance {
  // instance data
  string          name   = "";
  vector<frame3f> frames = {};
};

// Object.
struct object {
  // object data
  string    name     = "";
  frame3f   frame    = identity3x4f;
  shape*    shape    = nullptr;
  material* material = nullptr;
  instance* instance = nullptr;
  subdiv*   subdiv   = nullptr;
};

// Environment map.
struct environment {
  string   name         = "";
  frame3f  frame        = identity3x4f;
  vec3f    emission     = {0, 0, 0};
  texture* emission_tex = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct model {
  string               name         = "";
  vector<camera*>      cameras      = {};
  vector<object*>      objects      = {};
  vector<environment*> environments = {};
  vector<shape*>       shapes       = {};
  vector<subdiv*>      subdivs      = {};
  vector<texture*>     textures     = {};
  vector<material*>    materials    = {};
  vector<instance*>    instances    = {};
  ~model();
};

// add element to a scene
camera*      add_camera(model* scene, const string& name = "");
environment* add_environment(model* scene, const string& name = "");
object*      add_object(model* scene, const string& name = "");
instance*    add_instance(model* scene, const string& name = "");
material*    add_material(model* scene, const string& name = "");
shape*       add_shape(model* scene, const string& name = "");
subdiv*      add_subdiv(model* scene, const string& name = "");
texture*     add_texture(model* scene, const string& name = "");
object*      add_complete_object(model* scene, const string& name = "");

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto::sceneio {

// Progress callback called when loading.
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Load/save a scene in the supported formats. Throws on error.
// Calls the progress callback, if defined, as we process more data.
bool load_scene(const string& filename, model* scene, string& error,
    progress_callback progress_cb = {}, bool noparallel = false);
bool save_scene(const string& filename, const model* scene, string& error,
    progress_callback progress_cb = {}, bool noparallel = false);

// get named camera or default if name is empty
camera* get_camera(const model* scene, const string& name = "");

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Return scene statistics as list of strings.
vector<string> scene_stats(const model* scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(const model* scene, bool notextures = false);

// Return an approximate scene bounding box.
bbox3f compute_bounds(const model* scene);

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Apply subdivision and displacement rules.
void tesselate_subdivs(model* scene, progress_callback progress_cb = {});
void tesselate_subdiv(model* scene, subdiv* subdiv);

// Update node transforms. Eventually this will be deprecated as we do not
// support animation in this manner long term.
void update_transforms(
    model* scene, float time = 0, const string& anim_group = "");

// TODO: remove
inline vec3f eta_to_reflectivity(float eta) {
  return vec3f{((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1))};
}

}  // namespace yocto::sceneio

#endif
