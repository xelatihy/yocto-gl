//
// # Yocto/SceneIO: Tiny library for Yocto/Scene input and output
//
// Yocto/SceneIO provides loading and saving functionality for scenes
// in Yocto/GL. We support a simple to use YAML format, PLY, OBJ and glTF.
// The YAML serialization is a straight copy of the in-memory scene data.
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

#include "yocto_image.h"
#include "yocto_math.h"

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
  float   focus        = flt_max;
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
  // lobes
  vec3f emission     = {0, 0, 0};
  vec3f base         = {0, 0, 0};
  float specular     = 0;
  float roughness    = 0;
  float metallic     = 0;
  float ior          = 1.5;
  vec3f spectint     = {1, 1, 1};
  float coat         = 0;
  float transmission = 0;
  vec3f scattering   = {0, 0, 0};
  float phaseg       = 0;
  float radius       = 0.01;
  float opacity      = 1;
  bool  thin         = true;

  // textures
  int  emission_tex     = -1;
  int  base_tex         = -1;
  int  specular_tex     = -1;
  int  metallic_tex     = -1;
  int  roughness_tex    = -1;
  int  transmission_tex = -1;
  int  spectint_tex     = -1;
  int  scattering_tex   = -1;
  int  coat_tex         = -1;
  int  opacity_tex      = -1;
  int  normal_tex       = -1;
  bool gltf_textures    = false;  // glTF packed textures
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct sceneio_shape {
  // shape data
  string name      = "";
  frame3f frame = identity3x4f;
  vector<frame3f> instances = {};
  sceneio_material material = {};

  // instances

  // primitives
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};
};

// Subdiv data represented as indexed meshes of elements.
// May contain points, lines, triangles, quads or
// face-varying quads.
struct sceneio_subdiv {
  // shape data
  string name     = "";
  int shape = -1;

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

  // subdision properties
  int  subdivisions = 0;
  bool catmullclark = false;
  bool smooth       = false;

  // displacement information
  float displacement     = 0;
  int   displacement_tex = -1;
};

// Environment map.
struct sceneio_environment {
  string  name         = "";
  frame3f frame        = identity3x4f;
  vec3f   emission     = {0, 0, 0};
  int     emission_tex = -1;
};

// Node in a transform hierarchy.
struct sceneio_node {
  string        name        = "";
  int           parent      = -1;
  frame3f       local       = identity3x4f;
  vec3f         translation = {0, 0, 0};
  vec4f         rotation    = {0, 0, 0, 1};
  vec3f         scale       = {1, 1, 1};
  vector<float> weights     = {};
  int           camera      = -1;
  int           shape       = -1;
  int           instance    = -1;
  int           environment = -1;

  // compute properties
  vector<int> children = {};
};

// Keyframe data.
struct sceneio_animation {
  enum struct interpolation_type { linear, step, bezier };
  string                name          = "";
  string                filename      = "";
  string                group         = "";
  interpolation_type    interpolation = interpolation_type::linear;
  vector<float>         times         = {};
  vector<vec3f>         translations  = {};
  vector<vec4f>         rotations     = {};
  vector<vec3f>         scales        = {};
  vector<vector<float>> morphs        = {};
  vector<int>           targets       = {};
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct sceneio_model {
  string                      name         = "";
  vector<sceneio_camera>      cameras      = {};
  vector<sceneio_shape>       shapes       = {};
  vector<sceneio_subdiv>      subdivs      = {};
  vector<sceneio_texture>     textures     = {};
  vector<sceneio_environment> environments = {};
  vector<sceneio_node>        nodes        = {};
  vector<sceneio_animation>   animations   = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto {

// Load/save a scene in the supported formats. Throws on error.
void load_scene(
    const string& filename, sceneio_model& scene, bool noparallel = false);
void save_scene(const string& filename, const sceneio_model& scene,
    bool noparallel = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Return scene statistics as list of strings.
vector<string> scene_stats(const sceneio_model& scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
    const sceneio_model& scene, bool notextures = false);

// Return an approximate scene bounding box.
bbox3f compute_bounds(const sceneio_model& scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_subdiv(
    sceneio_model& scene, const sceneio_subdiv& subdiv, bool no_quads = false);

// Update node transforms. Eventually this will be deprecated as we do not
// support animation in this manner long term.
void update_transforms(
    sceneio_model& scene, float time = 0, const string& anim_group = "");

// TODO: remove
inline vec3f eta_to_reflectivity(float eta) {
  return vec3f{((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1))};
}

}  // namespace yocto

#endif
