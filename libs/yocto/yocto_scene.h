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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Handles to refer to scene elements
inline const int invalidid = -1;

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
  float   lens         = 0.050f;
  float   film         = 0.036f;
  float   aspect       = 1.500f;
  float   focus        = 10000;
  float   aperture     = 0;
};

// Texture data as array of float or byte pixels. Textures can be stored in
// linear or non linear color space.
struct scene_texture {
  int           width   = 0;
  int           height  = 0;
  bool          linear  = false;
  vector<vec4f> pixelsf = {};
  vector<vec4b> pixelsb = {};
};

// Material type
enum struct scene_material_type {
  // clang-format off
  matte, glossy, metallic, transparent, refractive, subsurface, volume, gltfpbr
  // clang-format on
};

// Enum labels
inline const auto scene_material_names = std::vector<std::string>{"matte",
    "glossy", "metallic", "transparent", "refractive", "subsurface", "volume",
    "gltfpbr"};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct scene_material {
  // material
  scene_material_type type         = scene_material_type::matte;
  vec3f               emission     = {0, 0, 0};
  vec3f               color        = {0, 0, 0};
  float               roughness    = 0;
  float               metallic     = 0;
  float               ior          = 1.5f;
  vec3f               scattering   = {0, 0, 0};
  float               scanisotropy = 0;
  float               trdepth      = 0.01f;
  float               opacity      = 1;

  // textures
  int emission_tex   = invalidid;
  int color_tex      = invalidid;
  int roughness_tex  = invalidid;
  int scattering_tex = invalidid;
  int normal_tex     = invalidid;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct scene_shape {
  // element data
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

// Shape data stored as a face-varying mesh
struct scene_fvshape {
  // element data
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
};

// Instance.
struct scene_instance {
  // instance data
  frame3f frame    = identity3x4f;
  int     shape    = invalidid;
  int     material = invalidid;
};

// Environment map.
struct scene_environment {
  // environment data
  frame3f frame        = identity3x4f;
  vec3f   emission     = {0, 0, 0};
  int     emission_tex = invalidid;
};

// Subdiv data represented as face-varying primitives where
// each vertex data has its own topology.
struct scene_subdiv {
  // face-varying primitives
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};

  // subdivision data
  int  subdivisions = 0;
  bool catmullclark = true;
  bool smooth       = true;

  // displacement data
  float displacement     = 0;
  int   displacement_tex = invalidid;

  // shape reference
  int shape = invalidid;
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
  vector<scene_camera>      cameras      = {};
  vector<scene_instance>    instances    = {};
  vector<scene_environment> environments = {};
  vector<scene_shape>       shapes       = {};
  vector<scene_texture>     textures     = {};
  vector<scene_material>    materials    = {};
  vector<scene_subdiv>      subdivs      = {};

  // names (this will be cleanup significantly later)
  vector<string> camera_names      = {};
  vector<string> texture_names     = {};
  vector<string> material_names    = {};
  vector<string> shape_names       = {};
  vector<string> instance_names    = {};
  vector<string> environment_names = {};
  vector<string> subdiv_names      = {};

  // copyright info preserve in IO
  string copyright = "";
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// CAMERA PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera.
ray3f eval_camera(
    const scene_camera& camera, const vec2f& image_uv, const vec2f& lens_uv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates a texture
vec4f eval_texture(const scene_texture& texture, const vec2f& uv,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);
vec4f eval_texture(const scene_model& scene, int texture, const vec2f& uv,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// pixel access
vec4f lookup_texture(
    const scene_texture& texture, int i, int j, bool as_linear = false);

// conversion from image
scene_texture image_to_texture(const color_image& image);

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATERIAL PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Material parameters evaluated at a point on the surface
struct material_point {
  scene_material_type type         = scene_material_type::gltfpbr;
  vec3f               emission     = {0, 0, 0};
  vec3f               color        = {0, 0, 0};
  float               opacity      = 1;
  float               roughness    = 0;
  float               metallic     = 0;
  float               ior          = 1;
  vec3f               density      = {0, 0, 0};
  vec3f               scattering   = {0, 0, 0};
  float               scanisotropy = 0;
  float               trdepth      = 0.01f;
};

// Eval material to obtain emission, brdf and opacity.
material_point eval_material(const scene_model& scene,
    const scene_material& material, const vec2f& texcoord,
    const vec4f& shape_color = {1, 1, 1, 1});

// check if a material is a delta
bool is_delta(const scene_material& material);
bool is_delta(const material_point& material);

// check if a material has a volume
bool is_volumetric(const scene_material& material);
bool is_volumetric(const material_point& material);
bool is_volumetric(const scene_model& scene, const scene_instance& instance);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const scene_shape& shape, int element, const vec2f& uv);
vec3f eval_normal(const scene_shape& shape, int element, const vec2f& uv);
vec3f eval_tangent(const scene_shape& shape, int element, const vec2f& uv);
vec2f eval_texcoord(const scene_shape& shape, int element, const vec2f& uv);
vec4f eval_color(const scene_shape& shape, int element, const vec2f& uv);
float eval_radius(const scene_shape& shape, int element, const vec2f& uv);

// Evaluate element normals
vec3f eval_element_normal(const scene_shape& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const scene_shape& shape);
void          compute_normals(vector<vec3f>& normals, const scene_shape& shape);

// An unevaluated location on a shape
struct shape_point {
  int   element = 0;
  vec2f uv      = {0, 0};
};

// Shape sampling
vector<float> sample_shape_cdf(const scene_shape& shape);
void          sample_shape_cdf(vector<float>& cdf, const scene_shape& shape);
shape_point   sample_shape(const scene_shape& shape, const vector<float>& cdf,
      float rn, const vec2f& ruv);
vector<shape_point> sample_shape(const scene_shape& shape,
    const vector<float>& cdf, int num_samples, uint64_t seed = 98729387);

// Conversions
scene_shape quads_to_triangles(const scene_shape& shape);
void        quads_to_triangles(scene_shape& result, const scene_shape& shape);

// Subdivision
scene_shape subdivide_shape(
    const scene_shape& shape, int subdivisions, bool catmullclark);

// Interpolate vertex data
vec3f eval_position(const scene_fvshape& shape, int element, const vec2f& uv);
vec3f eval_normal(const scene_fvshape& shape, int element, const vec2f& uv);
vec2f eval_texcoord(const scene_shape& shape, int element, const vec2f& uv);

// Evaluate element normals
vec3f eval_element_normal(const scene_fvshape& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const scene_fvshape& shape);
void compute_normals(vector<vec3f>& normals, const scene_fvshape& shape);

// Conversions
scene_shape fvshape_to_shape(
    const scene_fvshape& shape, bool as_triangles = false);
scene_fvshape shape_to_fvshape(const scene_shape& shape);

// Subdivision
scene_fvshape subdivide_fvshape(
    const scene_fvshape& shape, int subdivisions, bool catmullclark);

// Shape statistics
vector<string> shape_stats(const scene_shape& shape, bool verbose = false);
vector<string> fvshape_stats(const scene_fvshape& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// INSTANCE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate instance properties
vec3f eval_position(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec3f eval_element_normal(
    const scene_model& scene, const scene_instance& instance, int element);
vec3f eval_normal(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec2f eval_texcoord(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
    const scene_model& scene, const scene_instance& instance, int element);
vec3f eval_normalmap(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv);
vec3f eval_shading_normal(const scene_model& scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing);
vec4f eval_color(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv);

// Eval material to obtain emission, brdf and opacity.
material_point eval_material(const scene_model& scene,
    const scene_instance& instance, int element, const vec2f& uv);
// check if a material has a volume
bool is_volumetric(const scene_model& scene, const scene_instance& instance);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENVIRONMENT PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Environment
vec3f eval_environment(const scene_model& scene,
    const scene_environment& environment, const vec3f& direction);
vec3f eval_environment(const scene_model& scene, const vec3f& direction);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// compute scene bounds
bbox3f compute_bounds(const scene_model& scene);

// add missing elements
void add_camera(scene_model& scene);
void add_sky(scene_model& scene, float sun_angle = pif / 4);

// get named camera or default if name is empty
int find_camera(const scene_model& scene, const string& name);

// create a scene from a shape
scene_model make_shape_scene(const scene_shape& shape, bool add_sky = false);

// Return scene statistics as list of strings.
vector<string> scene_stats(const scene_model& scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
    const scene_model& scene, bool notextures = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

// Apply subdivision and displacement rules.
void tesselate_subdivs(scene_model& scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a plane.
scene_shape make_rect(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
scene_shape make_bulged_rect(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1},
    float radius = 0.3);
// Make a plane in the xz plane.
scene_shape make_recty(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
scene_shape make_bulged_recty(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1},
    float radius = 0.3);
// Make a box.
scene_shape make_box(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec3f& uvscale = {1, 1, 1});
scene_shape make_rounded_box(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec3f& uvscale = {1, 1, 1},
    float radius = 0.3);
// Make a quad stack
scene_shape make_rect_stack(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec2f& uvscale = {1, 1});
// Make a floor.
scene_shape make_floor(const vec2i& steps = {1, 1},
    const vec2f& scale = {10, 10}, const vec2f& uvscale = {10, 10});
scene_shape make_bent_floor(const vec2i& steps = {1, 1},
    const vec2f& scale = {10, 10}, const vec2f& uvscale = {10, 10},
    float bent = 0.5);
// Make a sphere.
scene_shape make_sphere(int steps = 32, float scale = 1, float uvscale = 1);
// Make a sphere.
scene_shape make_uvsphere(const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a sphere with slipped caps.
scene_shape make_capped_uvsphere(const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1}, float height = 0.3);
// Make a disk
scene_shape make_disk(int steps = 32, float scale = 1, float uvscale = 1);
// Make a bulged disk
scene_shape make_bulged_disk(
    int steps = 32, float scale = 1, float uvscale = 1, float height = 0.3);
// Make a uv disk
scene_shape make_uvdisk(const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv cylinder
scene_shape make_uvcylinder(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a rounded uv cylinder
scene_shape make_rounded_uvcylinder(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1},
    float radius = 0.3);

// Make a facevarying rect
scene_fvshape make_fvrect(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1});
// Make a facevarying box
scene_fvshape make_fvbox(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a facevarying sphere
scene_fvshape make_fvsphere(int steps = 32, float scale = 1, float uvscale = 1);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
scene_shape make_lines(const vec2i& steps = {4, 65536},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1},
    const vec2f& radius = {0.001f, 0.001f});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
scene_shape make_point(float radius = 0.001f);
scene_shape make_points(
    int num = 65536, float uvscale = 1, float radius = 0.001f);
scene_shape make_random_points(int num = 65536, const vec3f& size = {1, 1, 1},
    float uvscale = 1, float radius = 0.001f, uint64_t seed = 17);

// Predefined meshes
scene_shape   make_monkey(float scale = 1);
scene_shape   make_quad(float scale = 1);
scene_shape   make_quady(float scale = 1);
scene_shape   make_cube(float scale = 1);
scene_fvshape make_fvcube(float scale = 1);
scene_shape   make_geosphere(float scale = 1);

// Make a hair ball around a shape.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (strength/number)
// rotation: rotation added to hair (angle/strength)
scene_shape make_hair(const scene_shape& shape, const vec2i& steps = {8, 65536},
    const vec2f& length = {0.1f, 0.1f}, const vec2f& rad = {0.001f, 0.001f},
    const vec2f& noise = {0, 10}, const vec2f& clump = {0, 128},
    const vec2f& rotation = {0, 0}, int seed = 7);

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res and without texcoords and normals.
scene_shape points_to_spheres(
    const vector<vec3f>& vertices, int steps = 2, float scale = 0.01f);
scene_shape polyline_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);
scene_shape lines_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);

// Make a heightfield mesh.
scene_shape make_heightfield(const vec2i& size, const vector<float>& height);
scene_shape make_heightfield(const vec2i& size, const vector<vec4f>& color);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

// Make Cornell Box scene
void make_cornellbox(scene_model& scene);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARDS COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

using sceneio_scene       = scene_model;
using sceneio_camera      = scene_camera;
using sceneio_texture     = scene_texture;
using sceneio_material    = scene_material;
using sceneio_shape       = scene_shape;
using sceneio_instance    = scene_instance;
using sceneio_environment = scene_environment;

}  // namespace yocto

#endif
