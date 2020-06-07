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
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
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
  float displacement = 0;
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
  scene_texture* displacement_tex = nullptr;

  // [experimental] properties to drive subdiv and displacement
  int  subdivisions = 2;
  bool smooth       = true;
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
struct scene_subdiv {
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
struct scene_instance {
  // instance data
  string          name   = "";
  vector<frame3f> frames = {};
};

// Object.
struct scene_object {
  // object data
  string          name     = "";
  frame3f         frame    = identity3x4f;
  scene_shape*    shape    = nullptr;
  scene_material* material = nullptr;
  scene_instance* instance = nullptr;
  scene_subdiv*   subdiv   = nullptr;
};

// Environment map.
struct scene_environment {
  string         name         = "";
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {0, 0, 0};
  scene_texture* emission_tex = nullptr;
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
  vector<scene_subdiv*>      subdivs      = {};
  vector<scene_texture*>     textures     = {};
  vector<scene_material*>    materials    = {};
  vector<scene_instance*>    instances    = {};

  // additional information
  string name      = "";
  string copyright = "";

  // cleanup
  ~scene_model();
};

// add element to a scene
scene_camera*      add_camera(scene_model* scene, const string& name = "");
scene_environment* add_environment(scene_model* scene, const string& name = "");
scene_object*      add_object(scene_model* scene, const string& name = "");
scene_instance*    add_instance(scene_model* scene, const string& name = "");
scene_material*    add_material(scene_model* scene, const string& name = "");
scene_shape*       add_shape(scene_model* scene, const string& name = "");
scene_subdiv*      add_subdiv(scene_model* scene, const string& name = "");
scene_texture*     add_texture(scene_model* scene, const string& name = "");
scene_object* add_complete_object(scene_model* scene, const string& name = "");

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
bool load_scene(const string& filename, scene_model* scene, string& error,
    progress_callback progress_cb = {}, bool noparallel = false);
bool save_scene(const string& filename, const scene_model* scene, string& error,
    progress_callback progress_cb = {}, bool noparallel = false);

// get named camera or default if name is empty
scene_camera* get_camera(const scene_model* scene, const string& name = "");

// add a sky environment
void add_sky(scene_model* scene, float sun_angle = pif / 4);

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

// Return an approximate scene bounding box.
bbox3f compute_bounds(const scene_model* scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_subdivs(scene_model* scene, progress_callback progress_cb = {});
void tesselate_subdiv(scene_model* scene, scene_subdiv* subdiv);

}  // namespace yocto

#endif
