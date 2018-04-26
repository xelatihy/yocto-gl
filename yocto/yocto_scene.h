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
    float yfov = 2;                    // vertical field of view.
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float focus = 1;                   // focu distance
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
    std::vector<vec4b> ldr = {};  // ldr image
    std::vector<vec4f> hdr = {};  // hdr image
    void* gl_data = nullptr;      // unmanaged data for OpenGL viewer
};

// Texture information to use for lookup.
struct texture_info {
    texture* txt = nullptr;  // texture
    bool wrap_s = true;      // wrap in s coordinate
    bool wrap_t = true;      // wrop in t coordinate
    bool linear = true;      // linear interpolation
    bool mipmap = true;      // mipmapping
    float scale = 1;         // scale for occ, normal, bumps
};

// Material type.
enum struct material_type {
    specular_roughness = 0,   // microfacet model (OBJ)
    metallic_roughness = 1,   // base/metallic model (glTF)
    specular_glossiness = 2,  // sepcular/glossiness (glTF)
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on glTF for compatibility and adapted to OBJ.
// For lines, uses Kajija-Kay model. For points, a hacked up shading.
struct material {
    std::string name = "";      // name
    bool double_sided = false;  // double-sided rendering
    material_type type = material_type::specular_roughness;  // type

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f kd = {0, 0, 0};  // diffuse/base color
    vec3f ks = {0, 0, 0};  // specular color / metallic factor
    vec3f kr = {0, 0, 0};  // clear coat reflection
    vec3f kt = {0, 0, 0};  // transmission color
    float rs = 0.0001;     // roughness mapped as glTF
    float op = 1;          // opacity

    // textures
    texture_info ke_txt;    // emission texture
    texture_info kd_txt;    // diffuse texture
    texture_info ks_txt;    // specular texture
    texture_info kr_txt;    // clear coat reflection texture
    texture_info kt_txt;    // transmission texture
    texture_info rs_txt;    // roughness texture
    texture_info op_txt;    // opacity texture
    texture_info bump_txt;  // bump map texture (heighfield)
    texture_info disp_txt;  // displacement map texture (heighfield)
    texture_info norm_txt;  // normal texture
    texture_info occ_txt;   // occlusion texture
};

// Shape data represented as an indexed meshes of elements.
// May contain only element type (points/lines/triangles/quads/beziers).
struct shape {
    std::string name = "";  // name
    std::string path = "";  // path for glTF buffers

    // primitives
    std::vector<int> points;       // points
    std::vector<vec2i> lines;      // lines
    std::vector<vec3i> triangles;  // triangles
    std::vector<vec4i> quads;      // quads
    std::vector<vec4i> beziers;    // beziers

    // face-varying quad primitives
    std::vector<vec4i> quads_pos;       // pos indices
    std::vector<vec4i> quads_norm;      // norm indices
    std::vector<vec4i> quads_texcoord;  // texcoord indices

    // vertex data
    std::vector<vec3f> pos;        // positions
    std::vector<vec3f> norm;       // normals/tangents
    std::vector<vec2f> texcoord;   // texcoord coordinates
    std::vector<vec2f> texcoord1;  // second set of texture coordinates
    std::vector<vec4f> color;      // colors
    std::vector<float> radius;     // radia for lines/points
    std::vector<vec4f> tangsp;     // tangent space for triangles

    // tesselation data
    int subdivision = 0;        // subdivision [deprecated]
    bool catmullclark = false;  // catmull-clark [deprecated]

    // computed properties
    bbox3f bbox = invalid_bbox3f;      // boudning box
    std::vector<float> elem_cdf = {};  // element cdf for sampling
    bvh_tree* bvh = nullptr;           // bvh for ray intersection
    void* gl_data = nullptr;           // unmanaged data for OpenGL viewer

    // cleanup
    ~shape();
};

// Shape instance.
struct instance {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform frame
    shape* shp = nullptr;              // shape
    material* mat = nullptr;           // material

    // compute properties
    bbox3f bbox = invalid_bbox3f;  // boudning box
};

// Distance at which we set environment map positions.
const auto environment_distance = 1000000.0f;

// Envinonment map.
struct environment {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    vec3f ke = {0, 0, 0};              // emission color
    texture_info ke_txt = {};          // emission texture
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

// Light as either an instance or an environment.
struct light {
    instance* ist = nullptr;     // instance
    environment* env = nullptr;  // environment
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene {
    std::vector<camera*> cameras = {};            // cameras
    std::vector<shape*> shapes = {};              // shapes
    std::vector<instance*> instances = {};        // instances
    std::vector<material*> materials = {};        // materials
    std::vector<texture*> textures = {};          // textures
    std::vector<environment*> environments = {};  // environments

    std::vector<node*> nodes = {};            // node hierarchy [optional]
    std::vector<animation*> animations = {};  // animations [optional]

    // compute properties
    std::vector<light*> lights;
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
    bool preserve_quads = false, bool split_obj_shapes = true,
    bool skip_missing = true);
void save_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool preserve_obj_instances = false,
    bool gltf_separate_buffers = false, bool skip_missing = true);

// Loads/saves scene textures.
void load_textures(
    const std::string& filename, scene* scn, bool skip_missing = true);
void save_textures(
    const std::string& filename, const scene* scn, bool skip_missing = true);

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
void update_lights(scene* scn, bool do_shapes = true);
// Generate a distribution for sampling a shape uniformly based on area/length.
void update_shape_cdf(shape* shp);

// Updates/refits bvh.
void update_bvh(shape* shp, bool equalsize = false);
void update_bvh(scene* scn, bool do_shapes = true, bool equalsize = false);
void refit_bvh(shape* shp);
void refit_bvh(scene* scn, bool do_shapes = true);

// Subdivides shape elements.
void subdivide_shape_once(shape* shp, bool subdiv = false);
// Tesselate a shape into basic primitives.
void tesselate_shape(shape* shp, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);
// Tesselate scene shapes.
void tesselate_shapes(scene* scn, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);

// Add missing names, normals, tangents and hierarchy.
void add_missing_camera(scene* scn);
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
    int eid = 0;              // shape element index
    vec2f euv = zero2f;       // shape element coordinates
    float ray_t = 0;          // ray/point distance
};

// Intersects a ray with the scene.
scene_intersection intersect_ray(
    const scene* scn, const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv);
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv);
vec2f eval_texcoord(const shape* shp, int eid, const vec2f& euv);
vec4f eval_color(const shape* shp, int eid, const vec2f& euv);
float eval_radius(const shape* shp, int eid, const vec2f& euv);
vec4f eval_tangsp(const shape* shp, int eid, const vec2f& euv);
vec3f eval_tangsp(const shape* shp, int eid, const vec2f& euv, bool& left_handed);
// Shape element values.
vec3f eval_elem_norm(const shape* shp, int eid);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f eval_pos(const instance* ist, int eid, const vec2f& euv);
vec3f eval_norm(const instance* ist, int eid, const vec2f& euv);
vec2f eval_texcoord(const instance* ist, int eid, const vec2f& euv);
vec4f eval_color(const instance* ist, int eid, const vec2f& euv);
float eval_radius(const instance* ist, int eid, const vec2f& euv);
vec3f eval_tangsp(const instance* ist, int eid, const vec2f& euv, bool& left_handed);
// Instance element values.
vec3f eval_elem_norm(const instance* ist, int eid);

// Environment values interpolated using uv parametrization.
vec3f eval_pos(const environment* env, const vec2f& uv);
vec3f eval_norm(const environment* env, const vec2f& uv);
// Environment texture coordinates from uv parametrization.
vec2f eval_texcoord(const environment* env, const vec2f& uv);
// Evaluate uv parameters for an incoming direction.
vec2f eval_uv(const environment* env, const vec3f& w);

// Evaluate a texture.
vec4f eval_texture(const texture_info& info, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1});
inline vec4f eval_texture(const texture* txt, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1}) {
    auto info = texture_info();
    info.txt = (texture*)txt;
    return eval_texture(info, texcoord, srgb, def);
}
// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv);
// Synchronizes a camera aspect with image width and height. Set image
// values any one is 0 or less. Set camera aspect otherwise.
void sync_camera_aspect(const camera* cam, int& width, int& height);

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(const shape* shp,
    const std::vector<float>& cdf, float re, const vec2f& ruv);

// Sample an environment uniformly.
vec2f sample_environment(const environment* env, const vec2f& ruv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE ELEMENT CREATION
// -----------------------------------------------------------------------------
namespace ygl {

// make scene elements
camera* make_camera(const std::string& name, const vec3f& from, const vec3f& to,
    float yfov, float aspect);
texture* make_texture(const std::string& name, const std::string& path = "",
    int width = 0, int height = 0, const std::vector<vec4b>& ldr = {},
    const std::vector<vec4f>& hdr = {});
material* make_material(const std::string& name,
    const vec3f& kd = {0.2f, 0.2f, 0.2f}, const vec3f& ks = {0, 0, 0},
    float rs = 1);
instance* make_instance(const std::string& name, shape* shp = nullptr,
    material* mat = nullptr, const frame3f& frame = identity_frame3f);
node* make_node(const std::string& name, camera* cam = nullptr,
    instance* ist = nullptr, environment* env = nullptr,
    const frame3f& frame = identity_frame3f);
environment* make_environment(const std::string& name,
    const vec3f& ke = {1, 1, 1}, texture* ke_txt = nullptr,
    const frame3f& frame = identity_frame3f);

// make procedural scene elements
inline std::vector<std::string>& make_camera_types() {
    static auto names = std::vector<std::string>{"cam1", "cam2", "cam3"};
    return names;
}
camera* make_proc_camera(const std::string& name, const std::string& type);
inline std::vector<std::string>& proc_texture_types() {
    static auto names = std::vector<std::string>{"grid", "checker", "colored",
        "rcolored", "bump", "uv", "gamma", "noise", "ridge", "fbm",
        "turbulence", "grid_norm", "bump_norm", "gammaf", "sky"};
    return names;
}
texture* make_proc_texture(const std::string& name, const std::string& type,
    int resolution = 512, float scale = 8.0f, float sky_sunangle = pi / 4,
    float bump_scale = 4);
inline std::vector<std::string>& proc_material_types() {
    static auto names = std::vector<std::string>{"emission", "matte", "plastic",
        "metal", "glass", "transparent", "carpaint"};
    return names;
}
material* make_proc_material(const std::string& name, const std::string& type,
    const vec3f& color = {1, 1, 1}, float roughness = 1,
    const std::string& txt = "", const std::string& norm = "");
inline std::vector<std::string>& proc_shape_types() {
    static auto names = std::vector<std::string>{"floor", "quad", "cube",
        "cube_rounded", "sphere", "sphere_cube", "geodesic_sphere",
        "sphere_flipcap", "disk", "disk_quad", "disk_bulged", "cylinder",
        "cylinder_rounded", "cylindery", "cylindery_rounded", "suzanne",
        "cube_subdiv", "suzanne_subdiv", "fvcube_subdiv", "matball", "point",
        "pointscube", "hairball", "beziercircle"};
    return names;
}
shape* make_proc_shape(const std::string& name, const std::string& type,
    const vec3i& tesselation = zero3i, const vec3f& size = zero3f,
    const vec3f& uvsize = zero3f, float rounded = 0.75f, float radius = 0.001f);
instance* make_proc_instance(const std::string& name, const std::string& stype,
    const std::string& mtype, const frame3f& frame = identity_frame3f);

// Makes the Cornell Box scene.
scene* make_cornell_box_scene();
// Makes a simple scene with up to three objects lined up
scene* make_simple_scene(const std::vector<std::string>& shapes,
    const std::vector<std::string>& mats, const std::string& lights,
    bool nodes = false, const std::vector<std::string>& animations = {},
    const std::string& floor_mat = "matte_grid");
// Make a scene with random instances
scene* make_random_instances_scene(
    const vec2i& num, const bbox3f& bbox, uint64_t seed = 13);

}  // namespace ygl

#endif
