//
// # Yocto/Scene: Tiny C++ Library for 3D Scene storage and manipulation
//
// Yocto/Scene define a simple scene data structure useful to create quick demos
// for rendering applications. It is used for quick OpenGL previwing and to
// support path tracing.
//
// In Yocto/Scene, shapes are represented as indexed collections of points,
// lines, triangles, quads and bezier segments. Each shape may contain
// only one element type. Shapes are organized into a scene by creating shape
// instances, each its own transform. Materials are specified like in OBJ and
// glTF andinclude emission, base-metallic and diffuse-specular parametrization,
// normal, occlusion and displacement mapping. Finally, the scene containers
// cameras and environment maps. Quad support in shapes is experimental and
// mostly supported for loading and saving. Lights in Yocto/Scene are pointers
// to either instances or environments.
// The scene supports an optional node hierarchy with animation modeled on
// the glTF model.
//
// Shapes are parametrized as in Yocto/Bvh with (u,v) for triangles written
// w.r.t the (v1-v0) and (v2-v0) axis respetively. Quads are internally handled
// as pairs of two triangles v0,v1,v3 and v2,v3,v1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. This is equivalent to Intel's
// Embree.
//
// ## Usage
//
// 1. load a scene with `load_scene()` and save it with `save_scene()`,
//    currently supporting OBJ and glTF
// 2. add missing data with `add_XXX()` functions
// 3. use `compute_bbox()` to compute element bounds
// 4. can merge scene together with `merge_into()`
// 5. make scene elements with `make_XXX()` functions
// 6. make procedural elements and scenes with `make_proc_XXX()` functions
// 7. for ray-intersection and closest point queries, a BVH can be created with
//    `make_shape_bvh()`/`make_scene_bvh()` and udpated with `refit_shape_bvh()`
//    and `refit_scene_bvh()` -- see Yocto/Bvh for intersection and
//    closest point functions
// 8. compute interpolated values over scene elements with `eval_XXX()`
//    functions
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#ifndef _YGL_SCENE_H_
#define _YGL_SCENE_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

// OBJ support
#ifndef YGL_OBJ
#define YGL_OBJ 1
#endif

// glTF support
#ifndef YGL_GLTF
#define YGL_GLTF 1
#endif

#include "yocto_math.h"

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// forward declaration
struct bvh_tree;

// Camera.
struct camera {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    bool ortho = false;                // orthographic
    float width = 0.036f;              // film width (default to 35mm)
    float height = 0.024f;             // film height (default to 35mm)
    float focal = 0.050f;              // focal length (defaul to 50 mm)
    float focus = flt_max;             // focal distance (default to inf focus)
    float aperture = 0;                // lens aperture
    float near = 0.01f;                // near plane distance
    float far = 10000;                 // far plane distance
};

// Texture containing either an LDR or HDR image.
struct texture {
    std::string name = "";        // name
    std::string path = "";        // file path
    int width = 0;                // width
    int height = 0;               // height
    std::vector<vec4f> pxl = {};  // image pixels
    bool wrap_s = true;           // wrap in s coordinate
    bool wrap_t = true;           // wrop in t coordinate
    bool linear = true;           // linear interpolation
    bool mipmap = true;           // mipmapping
    float scale = 1;              // scale for occ, normal, bumps
    bool srgb = true;         // whether to use srgb when loading/saving to 8bit
    void* gl_data = nullptr;  // unmanaged data for OpenGL viewer
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on glTF for compatibility and adapted to OBJ.
// For lines, uses Kajija-Kay model. For points, a hacked up shading.
struct material {
    std::string name = "";       // name
    bool base_metallic = false;  // base-metallic parametrization (glTF)
    bool double_sided = false;   // double sided rendering

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f kd = {0, 0, 0};  // diffuse/base color
    vec3f ks = {0, 0, 0};  // specular color / metallic factor
    vec3f kt = {0, 0, 0};  // transmission color
    float rs = 0.0001;     // roughness mapped as glTF
    float op = 1;          // opacity
    bool fresnel = true;   // whether to use fresnel in reflections/transmission
    bool refract = false;  // whether to use use refraction in tranmission

    // textures
    texture* ke_txt;    // emission texture
    texture* kd_txt;    // diffuse texture
    texture* ks_txt;    // specular texture
    texture* kt_txt;    // transmission texture
    texture* rs_txt;    // roughness texture
    texture* op_txt;    // opacity texture
    texture* occ_txt;   // occlusion texture
    texture* bump_txt;  // bump map texture (heighfield)
    texture* disp_txt;  // displacement map texture (heighfield)
    texture* norm_txt;  // normal texture
};

// Shape data represented as an indexed meshes of elements.
// May contain either tringles, lines or a set of vertices.
struct shape {
    std::string name = "";  // name
    std::string path = "";  // path for glTF buffers

    // primitives
    std::vector<vec2i> lines;      // lines
    std::vector<vec3i> triangles;  // triangles

    // vertex data
    std::vector<vec3f> pos;       // positions
    std::vector<vec3f> norm;      // normals/tangents
    std::vector<vec2f> texcoord;  // texcoord coordinates
    std::vector<vec4f> color;     // colors
    std::vector<float> radius;    // radia for lines/points
    std::vector<vec4f> tangsp;    // tangent space for triangles

    // computed properties
    bbox3f bbox = invalid_bbox3f;      // boudning box
    std::vector<float> elem_cdf = {};  // element cdf for sampling
    bvh_tree* bvh = nullptr;           // bvh for ray intersection
    void* gl_data = nullptr;           // unmanaged data for OpenGL viewer

    // cleanup
    ~shape();
};

// Subdivision surface.
struct subdiv {
    std::string name = "";        // name
    int level = 0;                // subdivision level
    bool catmull_clark = true;    // catmull clark subdiv
    bool compute_normals = true;  // faceted subdivision

    // primitives
    std::vector<vec4i> quads_pos;       // quads for position
    std::vector<vec4i> quads_texcoord;  // quads for texture coordinates
    std::vector<vec4i> quads_color;     // quads for color

    // creases
    std::vector<vec3i> crease_pos;       // crease for position
    std::vector<vec3i> crease_texcoord;  // crease for texture coordinates

    // vertex data
    std::vector<vec3f> pos;       // positions
    std::vector<vec2f> texcoord;  // texcoord coordinates
    std::vector<vec4f> color;     // colors
};

// Shape instance.
struct instance {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform frame
    shape* shp = nullptr;              // shape
    material* mat = nullptr;           // material
    subdiv* sbd = nullptr;             // subdivision shape

    // compute properties
    bbox3f bbox = invalid_bbox3f;  // boudning box
};

// Envinonment map.
struct environment {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    vec3f ke = {0, 0, 0};              // emission color
    texture* ke_txt = {};              // emission texture

    // computed properties
    std::vector<float> elem_cdf;  // element cdf for sampling
};

// Node in a transform hierarchy.
struct node {
    std::string name = "";             // name
    node* parent = nullptr;            // parent
    frame3f frame = identity_frame3f;  // transform frame
    vec3f translation = {0, 0, 0};     // translation
    vec4f rotation = {0, 0, 0, 1};     // rotation
    vec3f scale = {1, 1, 1};           // scale
    std::vector<float> weights = {};   // morph weights
    camera* cam = nullptr;             // camera
    instance* ist = nullptr;           // instance
    environment* env = nullptr;        // environment

    // compute properties
    std::vector<node*> children = {};  // child nodes
};

// Keyframe type.
enum struct animation_type { linear, step, bezier };

// Keyframe data.
struct animation {
    std::string name;                              // name
    std::string path = "";                         // path for glTF buffer
    std::string group;                             // group
    animation_type type = animation_type::linear;  // type
    std::vector<float> times;                      // keyframe times
    std::vector<vec3f> translation;                // translation keyframes
    std::vector<vec4f> rotation;                   // rotation keyframes
    std::vector<vec3f> scale;                      // scale keyframes
    std::vector<std::vector<float>> weights;       // mprph weight keyframes
    std::vector<node*> targets;                    // target nodes
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene {
    std::string name;                             // name
    std::vector<camera*> cameras = {};            // cameras
    std::vector<shape*> shapes = {};              // shapes
    std::vector<subdiv*> subdivs = {};            // subdivs
    std::vector<instance*> instances = {};        // instances
    std::vector<material*> materials = {};        // materials
    std::vector<texture*> textures = {};          // textures
    std::vector<environment*> environments = {};  // environments

    std::vector<node*> nodes = {};            // node hierarchy [optional]
    std::vector<animation*> animations = {};  // animations [optional]

    // compute properties
    std::vector<instance*> lights;
    bbox3f bbox = invalid_bbox3f;  // boudning box
    bvh_tree* bvh = nullptr;

    // cleanup
    ~scene();
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE INPUT AND OUTPUT AND STATS
// -----------------------------------------------------------------------------
namespace ygl {

// Loads/saves a scene in OBJ and glTF formats.
scene* load_scene(const std::string& filename, bool load_textures = true,
    bool split_obj_shapes = true, bool skip_missing = true);
void save_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool preserve_obj_instances = false,
    bool preserve_obj_subdivs = true, bool gltf_separate_buffers = false,
    bool skip_missing = true);

// Print scene statistics.
void print_stats(scene* scn);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(scene* merge_into, scene* merge_from);

}  // namespace ygl

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Update the normals of a shape.
void update_normals(shape* shp);

// Update node transforms.
void update_transforms(
    scene* scn, float time = 0, const std::string& anim_group = "");
// Compute animation range.
vec2f compute_animation_range(
    const scene* scn, const std::string& anim_group = "");

// Computes shape/scene approximate bounds.
void udpate_bbox(shape* shp);
void update_bbox(scene* scn, bool do_shapes = true);

// Update lights.
void update_lights(
    scene* scn, bool do_shapes = true, bool do_environments = false);
// Generate a distribution for sampling a shape uniformly based on area/length.
void update_shape_cdf(shape* shp);
// Generate a distribution for sampling an environment texture uniformly
// based on angle and texture intensity.
void update_environment_cdf(environment* env);

// Updates/refits bvh.
void update_bvh(shape* shp, bool equalsize = false);
void update_bvh(scene* scn, bool do_shapes = true, bool equalsize = false);
void refit_bvh(shape* shp);
void refit_bvh(scene* scn, bool do_shapes = true);

// Updates tesselation.
void update_tesselation(const subdiv* sbd, shape* shp);
void update_tesselation(scene* scn);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(scene* scn);
void add_missing_normals(scene* scn);
void add_missing_tangent_space(scene* scn);
// Checks for validity of the scene.
std::vector<std::string> validate(const scene* scn, bool skip_textures = false);

}  // namespace ygl

// -----------------------------------------------------------------------------
// INTERSECTION, EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection.
struct scene_intersection {
    instance* ist = nullptr;  // instance or null for no intersection
    int ei = 0;               // shape element index
    vec2f uv = zero2f;        // shape element coordinates
    float dist = 0;           // ray/point distance
};

// Intersects a ray with the scene.
scene_intersection intersect_ray(
    const scene* scn, const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f eval_pos(const shape* shp, int ei, const vec2f& uv);
vec3f eval_norm(const shape* shp, int ei, const vec2f& uv);
vec2f eval_texcoord(const shape* shp, int ei, const vec2f& uv);
vec4f eval_color(const shape* shp, int ei, const vec2f& uv);
float eval_radius(const shape* shp, int ei, const vec2f& uv);
vec4f eval_tangsp(const shape* shp, int ei, const vec2f& uv);
vec3f eval_tangsp(const shape* shp, int ei, const vec2f& uv, bool& left_handed);
// Shape element values.
vec3f eval_elem_norm(const shape* shp, int ei);
vec4f eval_elem_tangsp(const shape* shp, int ei);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f eval_pos(const instance* ist, int ei, const vec2f& uv);
vec3f eval_norm(const instance* ist, int ei, const vec2f& uv);
vec2f eval_texcoord(const instance* ist, int ei, const vec2f& uv);
vec4f eval_color(const instance* ist, int ei, const vec2f& uv);
float eval_radius(const instance* ist, int ei, const vec2f& uv);
vec3f eval_tangsp(
    const instance* ist, int ei, const vec2f& uv, bool& left_handed);
// Instance element values.
vec3f eval_elem_norm(const instance* ist, int ei);
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const instance* ist, int ei, const vec2f& uv, const vec3f& o);

// Environment texture coordinates from the incoming direction.
vec2f eval_texcoord(const environment* env, const vec3f& i);
// Evaluate the incoming direction from the uv.
vec3f eval_direction(const environment* env, const vec2f& uv);
// Evaluate the environment emission.
vec3f eval_environment(const environment* env, const vec3f& i);

// Evaluate a texture.
vec4f eval_texture(const texture* info, const vec2f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const camera* cam);
float eval_camera_aspect(const camera* cam);
void set_camera_fovy(
    camera* cam, float fovy, float aspect, float width = 0.036f);

// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, const vec2i& imsize,
    const vec2f& puv, const vec2f& luv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f eval_emission(const instance* ist, int ei, const vec2f& uv);
vec3f eval_diffuse(const instance* ist, int ei, const vec2f& uv);
vec3f eval_specular(const instance* ist, int ei, const vec2f& uv);
vec3f eval_transmission(const instance* ist, int ei, const vec2f& uv);
float eval_roughness(const instance* ist, int ei, const vec2f& uv);
float eval_opacity(const instance* ist, int ei, const vec2f& uv);

// Material values packed into a convenience structure.
struct bsdf {
    vec3f kd = zero3f;     // diffuse
    vec3f ks = zero3f;     // specular
    vec3f kt = zero3f;     // transmission
    float rs = 1;          // roughness
    bool refract = false;  // whether to use refraction in transmission
};
bsdf eval_bsdf(const instance* ist, int ei, const vec2f& uv);
bool is_delta_bsdf(const bsdf& f);

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const shape* shp, float re, const vec2f& ruv);

// Sample an environment uniformly.
vec2f sample_environment(const environment* env, const vec2f& ruv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE ELEMENT CREATION
// -----------------------------------------------------------------------------
namespace ygl {

// make scene elements
camera* make_camera(const std::string& name, const frame3f& frame,
    float width = 0.036f, float height = 0.024f, float focal = 0.050f,
    float focus = flt_max, float aperture = 0);
camera* make_bbox_camera(const std::string& name, const bbox3f& bbox,
    float width = 0.036f, float height = 0.024f, float focal = 0.050f);
texture* make_texture(const std::string& name, const std::string& path = "",
    int width = 0, int height = 0, const std::vector<vec4f>& pixels = {},
    bool srgb = true);
material* make_material(const std::string& name,
    const vec3f& kd = {0.2f, 0.2f, 0.2f}, const vec3f& ks = {0, 0, 0},
    float rs = 1);
instance* make_instance(const std::string& name, shape* shp = nullptr,
    material* mat = nullptr, subdiv* sbd = nullptr,
    const frame3f& frame = identity_frame3f);
node* make_node(const std::string& name, camera* cam = nullptr,
    instance* ist = nullptr, environment* env = nullptr,
    const frame3f& frame = identity_frame3f);
environment* make_environment(const std::string& name,
    const vec3f& ke = {1, 1, 1}, texture* ke_txt = nullptr,
    const frame3f& frame = identity_frame3f);
animation* make_animation(const std::string& name,
    const std::vector<float>& times = {},
    const std::vector<vec3f>& translation = {},
    const std::vector<vec4f>& rotation = {},
    const std::vector<vec3f>& scale = {},
    const std::vector<node*>& targets = {}, bool bezier = false);

// make a scene
scene* make_scene(const std::string& name, const std::vector<camera*>& cams,
    const std::vector<instance*>& ists, const std::vector<environment*>& envs);
scene* make_scene(const std::string& name, const std::vector<camera*>& cams,
    const std::vector<instance*>& ists, const std::vector<environment*>& envs,
    const std::vector<node*>& ndes, const std::vector<animation*>& anms);

// example shapes
shape* make_floor_shape(
    const std::string& name, int tesselation = 0, float size = 40);
shape* make_quad_shape(
    const std::string& name, int tesselation = 0, float size = 2);
shape* make_sphere_shape(
    const std::string& name, int tesselation = 5, float size = 2);
shape* make_spherecube_shape(
    const std::string& name, int tesselation = 4, float size = 2);
shape* make_sphereflipcap_shape(
    const std::string& name, int tesselation = 5, float size = 2);
shape* make_cube_shape(
    const std::string& name, int tesselation = 0, float size = 2);
shape* make_cuberounded_shape(
    const std::string& name, int tesselation = 4, float size = 2);
shape* make_matball_shape(
    const std::string& name, int tesselation = 5, float size = 2);
shape* make_quadstack_shape(const std::string& name, int stack_tesselation = 4,
    int tesselation = 0, float size = 2);
shape* make_hairball_shape(const std::string& name, int hair_tesselation = 16,
    int tesselation = 2, float size = 2, const vec2f& len = {0.2f, 0.2f},
    const vec2f& noise = {0, 0}, const vec2f& clump = {0, 0},
    const vec2f& radius = {0.001f, 0.001f});

// example subdivs
subdiv* make_suzanne_subdiv(const std::string& name, int subdivision = 2);
subdiv* make_cube_subdiv(const std::string& name, int tesselation = 0,
    int subdivision = 4, float size = 2);
subdiv* make_quad_subdiv(const std::string& name, int tesselation = 0,
    int subdivision = 4, float size = 2);
subdiv* make_opencube_subdiv(const std::string& name, int tesselation = 0,
    int subdivision = 4, float size = 2);

// example materials
material* make_emission_material(const std::string& name, const vec3f& col,
    texture* txt = nullptr, texture* norm = nullptr);
material* make_matte_material(const std::string& name, const vec3f& col,
    texture* txt = nullptr, texture* norm = nullptr);
material* make_plastic_material(const std::string& name, const vec3f& col,
    float rs = 0.1f, texture* txt = nullptr, texture* norm = nullptr);
material* make_metal_material(const std::string& name, const vec3f& col,
    float rs = 0.1f, texture* txt = nullptr, texture* norm = nullptr);
material* make_glass_material(const std::string& name,
    const vec3f& col = {1, 1, 1}, float rs = 0, texture* txt = nullptr,
    texture* norm = nullptr);
material* make_solidglass_material(const std::string& name,
    const vec3f& col = {1, 1, 1}, float rs = 0, texture* txt = nullptr,
    texture* norm = nullptr);
material* make_transparent_material(const std::string& name, const vec3f& col,
    float op, texture* txt = nullptr, texture* norm = nullptr);

// example textures
texture* make_grid_texture(
    const std::string& name, int resolution = 512, int tile = 8);
texture* make_uvgrid_texture(
    const std::string& name, int resolution = 512, int tile = 8);
texture* make_bump_texture(
    const std::string& name, int resolution = 512, int tile = 16);
texture* make_bumpnorm_texture(const std::string& name, int resolution = 512,
    int tile = 16, float scale = 4.0f);
texture* make_sky_texture(
    const std::string& name, int resolution = 512, float skyangle = pi / 4);
texture* make_lights_texture(const std::string& name, int resolution = 512,
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pi / 4,
    float lwidth = pi / 16, float lheight = pi / 16);

// Makes the Cornell Box scene.
scene* make_cornellbox_scene(const std::string& name, bool envlight = false);
// Make a simple scene with three objects lined up.
scene* make_simple_scene(const std::string& name,
    const std::vector<instance*>& objs, bool envlight = false,
    const std::vector<animation*>& anms = {});
// Make a simple scene with single object on a floor.
scene* make_simple_scene(
    const std::string& name, shape* shp, material* mats, bool envlight = false);
// Make a simple scene with a single object.
scene* make_shape_scene(
    const std::string& name, shape* shp, material* mat, bool envlight = false);
// Make a simple scene with a single environment and a mirror ball.
scene* make_environment_scene(const std::string& name, environment* env);
// Make a scene with random instances
scene* make_random_instances_scene(const std::string& name, const vec2i& num,
    const bbox3f& bbox, uint64_t seed = 13);

}  // namespace ygl

#endif
