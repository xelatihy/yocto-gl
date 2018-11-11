//
// # Yocto/Scene: Tiny library for scene representation
//
//
// Yocto/Scene is a library to represent 3D scenes using a simple data-driven
// and value oriented design.
//
//
// # Simple scene representation
//
// Yocto/Scene define a simple scene data structure useful to create quick demos
// and as the repsetnation upon which the path tracer works.
//
// In Yocto scenes, shapes are represented as indexed collections of points,
// lines, triangles, quads and bezier segments. Each shape may contain
// only one element type. Shapes are organized into a scene by creating shape
// instances, each its own transform. Materials are specified like in OBJ and
// glTF and include emission, base-metallic and diffuse-specular
// parametrization, normal, occlusion and displacement mapping. Finally, the
// scene containers cameras and environment maps. Quad support in shapes is
// experimental and mostly supported for loading and saving. Lights in
// Yocto/Scene are pointers to either instances or environments. The scene
// supports an optional node hierarchy with animation modeled on the glTF model.
//
// 1. load a scene with Yocto/SceneIO,
// 2. use `compute_shape_box()/compute_scene_box()` to compute element bounds
// 3. can merge scene together with `merge_scene()`
// 4. for ray-intersection and closest point queries, a BVH can be created with
//    `build_shape_bvh()/build_scene_bvh()` and refit with
//    `refit_shape_bvh()/refit_scene_bvh()`
// 5. compute interpolated values over scene elements with `evaluate_XXX()`
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
//

#ifndef _YOCTO_SCENE_H_
#define _YOCTO_SCENE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Camera based on a simple lens model. The camera is placed using a frame.
// Camera projection is described in photorgaphics terms. In particular,
// we specify fil size (35mm by default), the focal lengthm the focus
// distance and the lens_aperture. All values are in meters.
struct yocto_camera {
    string  name           = "";
    frame3f frame          = identity_frame3f;
    bool    orthographic   = false;
    vec2f   film_size      = {0.036f, 0.024f};
    float   focal_length   = 0.050f;
    float   focus_distance = maxf;
    float   lens_aperture  = 0;
};

// Texture containing either an LDR or HDR image. Textures are rendered
// using linear interpolation (unless `no_interoilation` is set) and
// weith tiling (unless `clamp_to_edge` is set). HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB. The latter
// conversion can be disabled with `ldr_as_linear` for example to render
// normal maps.
struct yocto_texture {
    string       name             = "";
    string       filename         = "";
    image<vec4f> hdr_image        = {};
    image<vec4b> ldr_image        = {};
    bool         clamp_to_edge    = false;
    bool         no_interpolation = false;
    float        height_scale     = 1;
    bool         ldr_as_linear    = false;
    bool         has_opacity      = false;
};

// Volumetric texture containing a float only volume data. See texture
// above for other propoerties.
struct yocto_voltexture {
    string        name             = "";
    string        filename         = "";
    volume<float> volume_data      = {};
    bool          clamp_to_edge    = false;
    bool          no_interpolation = false;
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct yocto_material {
    string name          = "";
    bool   base_metallic = false;  // base-metallic parametrization
    bool   gltf_textures = false;  // glTF packed textures
    bool   double_sided  = false;  // double sided rendering

    // base values
    vec3f emission     = {0, 0, 0};
    vec3f diffuse      = {0, 0, 0};
    vec3f specular     = {0, 0, 0};
    vec3f transmission = {0, 0, 0};
    float roughness    = 0.0001;
    float opacity      = 1;
    bool  fresnel      = true;
    bool  refract      = false;

    // textures
    int emission_texture     = -1;
    int diffuse_texture      = -1;
    int specular_texture     = -1;
    int transmission_texture = -1;
    int roughness_texture    = -1;
    int opacity_texture      = -1;
    int occlusion_texture    = -1;
    int bump_texture         = -1;
    int displacement_texture = -1;
    int normal_texture       = -1;

    // volume properties
    // albedo = scattering / (absorption + scattering)
    // density = absorption + scattering
    vec3f volume_emission = {0, 0, 0};
    vec3f volume_albedo   = {0, 0, 0};
    vec3f volume_density  = {0, 0, 0};
    float volume_phaseg   = 0;

    // volume textures
    int volume_density_texture = -1;
};

// Shape data represented as an indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct yocto_shape {
    // shape data
    string name     = "";
    string filename = "";
    int    material = -1;

    // subdision properties
    int  subdivision_level      = 0;
    bool catmull_clark          = false;
    bool compute_vertex_normals = false;

    // primitives
    vector<int>   points    = {};
    vector<vec2i> lines     = {};
    vector<vec3i> triangles = {};
    vector<vec4i> quads     = {};

    // vertex data
    vector<vec3f> positions     = {};
    vector<vec3f> normals       = {};
    vector<vec2f> texturecoords = {};
    vector<vec4f> colors        = {};
    vector<float> radius        = {};
    vector<vec4f> tangentspaces = {};
};

// Shape data represented as an indexed meshes of elements with face-varying
// data. Each face maintains different topologies for positions, normals and
// texture coordinates.
struct yocto_surface {
    // shape data
    string      name      = "";
    string      filename  = "";
    vector<int> materials = {};

    // subdision properties
    int  subdivision_level      = 0;
    bool catmull_clark          = false;
    bool compute_vertex_normals = false;

    // face-varying primitives
    vector<vec4i> quads_positions     = {};
    vector<vec4i> quads_normals       = {};
    vector<vec4i> quads_texturecoords = {};
    vector<int>   quads_materials     = {};

    // vertex data
    vector<vec3f> positions     = {};
    vector<vec3f> normals       = {};
    vector<vec2f> texturecoords = {};
};

// Instance of a visible object in the scene. For now, this can be either
// a shape or a surface.
struct yocto_instance {
    string  name    = "";
    frame3f frame   = identity_frame3f;
    int     shape   = -1;
    int     surface = -1;
};

// Environment map.
struct yocto_environment {
    string  name             = "";
    frame3f frame            = identity_frame3f;
    vec3f   emission         = {0, 0, 0};
    int     emission_texture = -1;
};

// Node in a transform hierarchy.
struct yocto_scene_node {
    string        name        = "";
    int           parent      = -1;
    frame3f       local       = identity_frame3f;
    vec3f         translation = {0, 0, 0};
    vec4f         rotation    = {0, 0, 0, 1};
    vec3f         scale       = {1, 1, 1};
    vector<float> weights     = {};
    int           camera      = -1;
    int           instance    = -1;
    int           environment = -1;

    // compute properties
    vector<int> children = {};
};

// Keyframe type.
enum struct yocto_interpolation_type { linear, step, bezier };

// Keyframe data.
struct yocto_animation {
    string                   name               = "";
    string                   filename           = "";
    string                   animation_group    = "";
    yocto_interpolation_type interpolation_type = yocto_interpolation_type::linear;
    vector<float>            keyframes_times    = {};
    vector<vec3f>            translation_keyframes   = {};
    vector<vec4f>            rotation_keyframes      = {};
    vector<vec3f>            scale_keyframes         = {};
    vector<vector<float>>    morph_weights_keyframes = {};
    vector<int>              node_targets            = {};
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct yocto_scene {
    string                    name         = "";
    vector<yocto_camera>      cameras      = {};
    vector<yocto_shape>       shapes       = {};
    vector<yocto_surface>     surfaces     = {};
    vector<yocto_instance>    instances    = {};
    vector<yocto_material>    materials    = {};
    vector<yocto_texture>     textures     = {};
    vector<yocto_environment> environments = {};
    vector<yocto_voltexture>  voltextures  = {};
    vector<yocto_scene_node>  nodes        = {};
    vector<yocto_animation>   animations   = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print scene statistics.
void print_stats(const yocto_scene& scene);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_scene(yocto_scene& merge_into, const yocto_scene& merge_from);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Update node transforms.
void update_transforms(
    yocto_scene& scene, float time = 0, const string& anim_group = "");

// Compute animation range.
vec2f compute_animation_range(
    const yocto_scene& scene, const string& anim_group = "");

// Computes shape/scene approximate bounds.
bbox3f compute_shape_bounds(const yocto_shape& shape);
bbox3f compute_scene_bounds(const yocto_scene& scene);

// Compute shape vertex normals
vector<vec3f> compute_shape_normals(const yocto_shape& shape);
vector<vec3f> compute_surface_normals(const yocto_surface& surface);

// Updates/refits bvh.
bvh_shape make_shape_bvh(
    const yocto_shape& shape, const build_bvh_options& options = {});
bvh_shape make_surface_bvh(
    const yocto_surface& surface, const build_bvh_options& options = {});
bvh_scene make_scene_bvh(
    const yocto_scene& scene, const build_bvh_options& options = {});
void refit_shape_bvh(const yocto_shape& shape, bvh_shape& bvh);
void refit_surface_bvh(const yocto_surface& surface, bvh_shape& bvh);
void refit_scene_bvh(const yocto_scene& scene, bvh_scene& bvh);

// Apply subdivision and displacement rules.
void tesselate_shapes_and_surfaces(yocto_scene& scene);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(yocto_scene& scene);
void add_missing_normals(yocto_scene& scene);
void add_missing_tangent_space(yocto_scene& scene);
void add_missing_materials(yocto_scene& scene);
void add_missing_cameras(yocto_scene& scene);

// Add a sky environment
void add_sky_environment(yocto_scene& scene, float sun_angle = pif / 4);

// Checks for validity of the scene.
void log_validation_errors(const yocto_scene& scene, bool skip_textures = false);

// Queries on objects
bool is_shape_face_varying(const yocto_shape& shape);

// Scene intersection. Upron intersection we set the instance pointer,
// the shape element_id and element_uv and the inetrsection distance.
struct scene_intersection {
    int   instance_id = -1;
    int   element_id  = -1;
    vec2f element_uv  = zero2f;
    float distance    = maxf;
};

// Intersects a ray with the scene.
scene_intersection intersect_scene(const yocto_scene& scene,
    const bvh_scene& bvh, const ray3f& ray, bool find_any = false);
// Intersects a ray with a scene instance.
scene_intersection intersect_scene(const yocto_scene& scene, int instance_id,
    const bvh_scene& bvh, const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f evaluate_shape_position(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec3f evaluate_shape_normal(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec2f evaluate_shape_texturecoord(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec4f evaluate_shape_color(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
float evaluate_shape_radius(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec4f evaluate_shape_tangentspace(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec3f evaluate_shape_tangentspace(const yocto_shape& shape, int element_id,
    const vec2f& element_uv, bool& left_handed);
// Shape element values.
vec3f evaluate_shape_element_normal(const yocto_shape& shape, int element_id);
vec4f evaluate_shape_element_tangentspace(
    const yocto_shape& shape, int element_id);

// Sample a shape element based on area/length.
vector<float>    compute_shape_elements_cdf(const yocto_shape& shape);
pair<int, vec2f> sample_shape_element(const yocto_shape& shape,
    const vector<float>& elem_cdf, float re, const vec2f& ruv);
float            sample_shape_element_pdf(const yocto_shape& shape,
               const vector<float>& elem_cdf, int element_id, const vec2f& element_uv);

// Surface values interpolated using barycentric coordinates.
vec3f evaluate_surface_position(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
vec3f evaluate_surface_normal(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
vec2f evaluate_surface_texturecoord(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
// Surface element values.
vec3f evaluate_surface_element_normal(const yocto_surface& shape, int element_id);
// Per-element material.
int get_surface_element_material(const yocto_surface& surface, int element_id);

// Sample a surface element based on area.
vector<float>    compute_surface_elements_cdf(const yocto_surface& surface);
pair<int, vec2f> sample_surface_element(const yocto_surface& surface,
    const vector<float>& elem_cdf, float re, const vec2f& ruv);
float            sample_surface_element_pdf(const yocto_surface& surface,
               const vector<float>& elem_cdf, int element_id, const vec2f& element_uv);

// Evaluate a texture.
vec2i evaluate_texture_size(const yocto_texture& texture);
vec4f lookup_texture(const yocto_texture& texture, const vec2i& ij);
vec4f evaluate_texture(const yocto_texture& texture, const vec2f& texcoord);
float lookup_voltexture(const yocto_voltexture& texture, const vec3i& ijk);
float evaluate_voltexture(const yocto_voltexture& texture, const vec3f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float get_camera_fovx(const yocto_camera& camera);
float get_camera_fovy(const yocto_camera& camera);
float get_camera_aspect(const yocto_camera& camera);
vec2i get_camera_image_size(const yocto_camera& camera, const vec2i& size);
void  set_camera_fovy(
     yocto_camera& camera, float fovy, float aspect, float width = 0.036f);
// Sets camera field of view to enclose all the bbox. Camera view direction
// fiom size and forcal lemgth can be overridden if we pass non zero values.
void set_camera_view(yocto_camera& camera, const bbox3f& bbox,
    const vec3f& view_direction = zero3f, const vec2f& film = zero2f,
    float focal = 0);

// Generates a ray from a camera image coordinate and lens coordinates.
ray3f evaluate_camera_ray(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv);
// Generates a ray from a camera for pixel `image_ij`, the image size,
// the sub-pixel coordinates `pixel_uv` and the lens coordinates `lens_uv`
// and the image resolution `image_size`.
ray3f evaluate_camera_ray(const yocto_camera& camera, const vec2i& image_ij,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv);
// Generates a ray from a camera for pixel index `idx`, the image size,
// the sub-pixel coordinates `pixel_uv` and the lens coordinates `lens_uv`.
ray3f evaluate_camera_ray(const yocto_camera& camera, int idx,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f evaluate_material_emission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_diffuse(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_specular(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_transmission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
float evaluate_material_roughness(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
float evaluate_material_opacity(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
// Query material properties
bool is_material_emissive(const yocto_material& material);

// Material values packed into a convenience structure.
struct microfacet_brdf {
    vec3f diffuse      = zero3f;
    vec3f specular     = zero3f;
    vec3f transmission = zero3f;
    float roughness    = 1;
    bool  refract      = false;
};
microfacet_brdf evaluate_material_brdf(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
bool            is_brdf_delta(const microfacet_brdf& f);
bool            is_brdf_zero(const microfacet_brdf& f);

// Check volume properties.
bool is_material_volume_homogeneus(const yocto_material& vol);
bool is_material_volume_colored(const yocto_material& vol);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f evaluate_instance_position(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec2f evaluate_instance_texturecoord(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec4f evaluate_instance_color(
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_tangentspace(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    bool& left_handed);
// Instance element values.
vec3f evaluate_instance_element_normal(
    const yocto_scene& scene, const yocto_instance& instance, int element_id);
// Shading normals including material perturbations.
vec3f evaluate_instance_shading_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    const vec3f& o);

// Material values
int   get_instance_material_id(const yocto_scene& scene,
      const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_emission(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
float evaluate_instance_opacity(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
bool  is_instance_emissive(
     const yocto_scene& scene, const yocto_instance& instance);

// <aterial brdf
microfacet_brdf evaluate_instance_brdf(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);

// Environment texture coordinates from the incoming direction.
vec2f evaluate_environment_texturecoord(
    const yocto_environment& environment, const vec3f& direction);
// Evaluate the incoming direction from the element_uv.
vec3f evaluate_environment_direction(
    const yocto_environment& environment, const vec2f& environment_uv);
// Evaluate the environment emission.
vec3f evaluate_environment_emission(const yocto_scene& scene,
    const yocto_environment& environment, const vec3f& direction);
// Evaluate all environment emission.
vec3f evaluate_environment_emission(
    const yocto_scene& scene, const vec3f& direction);

// Sample an environment based on either texel values of uniform
vector<float> compute_environment_texels_cdf(
    const yocto_scene& scene, const yocto_environment& environment);
vec3f sample_environment_direction(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    float re, const vec2f& ruv);
float sample_environment_direction_pdf(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    const vec3f& direction);

}  // namespace yocto

#endif
