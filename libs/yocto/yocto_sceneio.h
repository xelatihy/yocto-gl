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
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Namespace aliases
namespace scn = yocto::sceneio;
namespace img = yocto::image;

// Math defitions
using math::bbox3f;
using math::byte;
using math::frame3f;
using math::identity3x3f;
using math::identity3x4f;
using math::mat4f;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4b;
using math::vec4f;
using math::vec4i;
using math::zero2f;
using math::zero2i;
using math::zero3f;
using math::zero3i;

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

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
  std::string name         = "";
  frame3f     frame        = identity3x4f;
  bool        orthographic = false;
  float       lens         = 0.050;
  float       film         = 0.036;
  float       aspect       = 1.500;
  float       focus        = 10000;
  float       aperture     = 0;
};

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
struct texture {
  std::string       name    = "";
  img::image<vec3f> colorf  = {};
  img::image<vec3b> colorb  = {};
  img::image<float> scalarf = {};
  img::image<byte>  scalarb = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct material {
  // material data
  std::string name = "";

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
  scn::texture* emission_tex     = nullptr;
  scn::texture* color_tex        = nullptr;
  scn::texture* specular_tex     = nullptr;
  scn::texture* metallic_tex     = nullptr;
  scn::texture* roughness_tex    = nullptr;
  scn::texture* transmission_tex = nullptr;
  scn::texture* spectint_tex     = nullptr;
  scn::texture* scattering_tex   = nullptr;
  scn::texture* coat_tex         = nullptr;
  scn::texture* opacity_tex      = nullptr;
  scn::texture* normal_tex       = nullptr;
  scn::texture* displacement_tex = nullptr;

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
  std::string name = "";

  // primitives
  std::vector<int>   points    = {};
  std::vector<vec2i> lines     = {};
  std::vector<vec3i> triangles = {};
  std::vector<vec4i> quads     = {};

  // vertex data
  std::vector<vec3f> positions = {};
  std::vector<vec3f> normals   = {};
  std::vector<vec2f> texcoords = {};
  std::vector<vec3f> colors    = {};
  std::vector<float> radius    = {};
  std::vector<vec4f> tangents  = {};
};

// Subdiv data represented as indexed meshes of elements.
// May contain points, lines, triangles, quads or
// face-varying quads.
struct subdiv {
  // shape data
  std::string name = "";

  // face-varying primitives
  std::vector<vec4i> quadspos      = {};
  std::vector<vec4i> quadsnorm     = {};
  std::vector<vec4i> quadstexcoord = {};

  // vertex data
  std::vector<vec3f> positions = {};
  std::vector<vec3f> normals   = {};
  std::vector<vec2f> texcoords = {};
};

// Instance data.
struct instance {
  // instance data
  std::string          name   = "";
  std::vector<frame3f> frames = {};
};

// Object.
struct object {
  // object data
  std::string    name     = "";
  frame3f        frame    = identity3x4f;
  scn::shape*    shape    = nullptr;
  scn::material* material = nullptr;
  scn::instance* instance = nullptr;
  scn::subdiv*   subdiv   = nullptr;
};

// Environment map.
struct environment {
  std::string   name         = "";
  frame3f       frame        = identity3x4f;
  vec3f         emission     = {0, 0, 0};
  scn::texture* emission_tex = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct model {
  // scene elements
  std::vector<scn::camera*>      cameras      = {};
  std::vector<scn::object*>      objects      = {};
  std::vector<scn::environment*> environments = {};
  std::vector<scn::shape*>       shapes       = {};
  std::vector<scn::subdiv*>      subdivs      = {};
  std::vector<scn::texture*>     textures     = {};
  std::vector<scn::material*>    materials    = {};
  std::vector<scn::instance*>    instances    = {};

  // additional information
  std::string name      = "";
  std::string copyright = "";

  // cleanup
  ~model();
};

// add element to a scene
scn::camera*      add_camera(scn::model* scene, const std::string& name = "");
scn::environment* add_environment(
    scn::model* scene, const std::string& name = "");
scn::object*   add_object(scn::model* scene, const std::string& name = "");
scn::instance* add_instance(scn::model* scene, const std::string& name = "");
scn::material* add_material(scn::model* scene, const std::string& name = "");
scn::shape*    add_shape(scn::model* scene, const std::string& name = "");
scn::subdiv*   add_subdiv(scn::model* scene, const std::string& name = "");
scn::texture*  add_texture(scn::model* scene, const std::string& name = "");
scn::object*   add_complete_object(
      scn::model* scene, const std::string& name = "");

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto::sceneio {

// Progress callback called when loading.
using progress_callback =
    std::function<void(const std::string& message, int current, int total)>;

// Load/save a scene in the supported formats. Throws on error.
// Calls the progress callback, if defined, as we process more data.
bool load_scene(const std::string& filename, scn::model* scene,
    std::string& error, progress_callback progress_cb = {},
    bool noparallel = false);
bool save_scene(const std::string& filename, const scn::model* scene,
    std::string& error, progress_callback progress_cb = {},
    bool noparallel = false);

// get named camera or default if name is empty
scn::camera* get_camera(const scn::model* scene, const std::string& name = "");

// add a sky environment
void add_sky(scn::model* scene, float sun_angle = math::pif / 4);

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Make Cornell Box scene
void make_cornellbox(scn::model* scene);

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Return scene statistics as list of strings.
std::vector<std::string> scene_stats(
    const scn::model* scene, bool verbose = false);
// Return validation errors as list of strings.
std::vector<std::string> scene_validation(
    const scn::model* scene, bool notextures = false);

// Return an approximate scene bounding box.
bbox3f compute_bounds(const scn::model* scene);

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Apply subdivision and displacement rules.
void tesselate_subdivs(scn::model* scene, progress_callback progress_cb = {});
void tesselate_subdiv(scn::model* scene, scn::subdiv* subdiv);

}  // namespace yocto::sceneio

#endif
