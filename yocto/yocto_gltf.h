///
/// # Yocto/glTF: High-level glTF support for Yocto/GL
///
/// Extension of Yocto/GL that provides a complete high-level interface for
/// Khronos glTF format. Uses the low-level parser from Yocto/GL.
/// The high-level interface that is more useful if an aplication needs direct
/// access to shapes, textures and animations and does not want to deal with
/// all the intricacies of the format.
///
/// Known limitations are: (1) skinning matrices are always in world space
/// (waiting for the spec to be updated); (2) spline based animation is not
/// fully implemented yet (waiting for official demo models).
///
/// This library depends on `yocto_gl.h` and its dependencies.
///
///
/// ## Usage
///
/// 1. load a group of scens with `load_scenes()`
/// 2. look at the `scene` data structures for access to individual elements
/// 3. to support animation, use `update_animated_transforms()`
/// 4. to support skinning, use `get_skin_transforms()`
/// 5. for morphing, use `compute_morphing_deformation()`
/// 6. can also manipulate the scene by adding missing data with `add_XXX()`
///    functions
/// 7. for rendering scenes, use `get_scene_cameras()` and
///    `get_scene_instances()` that avoid the need for explicirtly walking
///    the glTF node hierarchy
/// 8. use `save_scenes()` to write the data to disk
/// 9. use `convert_to_specgloss()` to convert materials to spec-gloss
///
///
/// ## History
///
/// - v 0.1.0: initial release after refactoring
///

// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YGLTF_H_
#define _YGLTF_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_gl.h"

// -----------------------------------------------------------------------------
// HIGH-LEVEL INTERFACE FOR KHRONOS GLTF
// -----------------------------------------------------------------------------
namespace ygl {

/// Camera
struct gltf_camera {
    /// name
    std::string name = "";
    /// orthographic
    bool ortho = false;
    /// aspect ratio
    float aspect = 1;
    /// vertical fov (perspective) or size (orthographic)
    float yfov = 2 * std::atan(0.5f);
    /// near plane (0 for default)
    float near = 0;
    /// far plane (0 for default)
    float far = 0;
    /// focus distance (extension not implemented yet)
    float focus = 1;
    /// lens aperture (extension not implemented yet)
    float aperture = 0;
};

/// Texture
struct gltf_texture {
    /// name
    std::string name = "";
    /// path
    std::string path = "";
    /// 8-bit data
    image4b ldr;
    /// float data
    image4f hdr;

    /// get texture width
    int width() const {
        if (!ldr.empty()) return ldr.width();
        if (!hdr.empty()) return hdr.width();
        return 0;
    }
    /// get texture height
    int height() const {
        if (!ldr.empty()) return ldr.height();
        if (!hdr.empty()) return hdr.height();
        return 0;
    }
};

/// Texture wrap mode
enum struct gltf_texture_wrap {
    /// repeat
    repeat = 1,
    /// clamp
    clamp = 2,
    /// mirror
    mirror = 3,
};

/// Texture filter mode
enum struct gltf_texture_filter {
    /// linear
    linear = 1,
    /// nearest
    nearest = 2,
    /// linear mipmap linear
    linear_mipmap_linear = 3,
    /// nearest mipmap nearest
    nearest_mipmap_nearest = 4,
    /// linear mipmap nearest
    linear_mipmap_nearest = 5,
    /// nearest mipmap linear
    nearest_mipmap_linear = 6,
};

/// Texture information
struct gltf_texture_info {
    /// wrap mode for s coordinate
    gltf_texture_wrap wrap_s = gltf_texture_wrap::repeat;
    /// wrap mdoe for t coordinate
    gltf_texture_wrap wrap_t = gltf_texture_wrap::repeat;
    /// filter mode
    gltf_texture_filter filter_mag = gltf_texture_filter::linear;
    /// filter mode
    gltf_texture_filter filter_min = gltf_texture_filter::linear_mipmap_linear;

    /// texture strength (occlusion and normal)
    float scale = 1;

    /// check if it is default
    bool is_default() const {
        return wrap_s == gltf_texture_wrap::repeat &&
               wrap_t == gltf_texture_wrap::repeat &&
               filter_mag == gltf_texture_filter::linear &&
               filter_min == gltf_texture_filter::linear_mipmap_linear;
    }
};

/// Material PBR metallic roughness
struct gltf_material_metallic_roughness {
    /// base color
    vec3f base = {0, 0, 0};
    /// opacity
    float opacity = 1;
    /// metallic factor
    float metallic = 0;
    /// metallic roughness
    float roughness = 0;
    /// base texture (kb.x, kb.y, kb.z, op)
    gltf_texture* base_txt = nullptr;
    /// metallic-roughness texture (n/a, roughness, metallic, n/a)
    gltf_texture* metallic_txt = nullptr;
    /// texture information for base_txt
    gltf_texture_info* base_txt_info = nullptr;
    /// texture information for metallic_txt
    gltf_texture_info* metallic_txt_info = nullptr;

    // cleanup
    ~gltf_material_metallic_roughness() {
        if (base_txt_info) delete base_txt_info;
        if (metallic_txt_info) delete metallic_txt_info;
    }
};

/// Material PBR specular glossiness
struct gltf_material_specular_glossiness {
    /// diffuse color and opacity
    vec3f diffuse = {0, 0, 0};
    /// opacity
    float opacity = 1;
    /// specular color (spec.x, spec.y, spec.z, opacity)
    vec3f specular = {0, 0, 0};
    /// specular glossiness
    float glossiness = 1;
    /// diffuse texture (diff.x, diff.y, diff.z, opacity)
    gltf_texture* diffuse_txt = nullptr;
    /// specular-glossiness texture (spec.x, spec.y, spec.z, gloss)
    gltf_texture* specular_txt = nullptr;
    /// texture information for base_txt
    gltf_texture_info* diffuse_txt_info = nullptr;
    /// texture information for metallic_txt
    gltf_texture_info* specular_txt_info = nullptr;

    // cleanup
    gltf_material_specular_glossiness() {
        if (diffuse_txt_info) delete diffuse_txt_info;
        if (specular_txt_info) delete specular_txt_info;
    }
};

/// Material
///
/// glTF 2.0 has two physically-based aterial models: pbrMetallicRoughness
/// and pbrSpecularGlossiness, the latter as an extension. Here we support both
/// by including which one is defined. While it would have been more appropriate
/// to convert them, this requires a rewrite of texture data which w prefer to
/// avoid.
struct gltf_material {
    /// name
    std::string name = "";

    // common properties --------------------
    /// emission color
    vec3f emission = {0, 0, 0};
    /// emissive texture reference
    gltf_texture* emission_txt = nullptr;
    /// texture information for normal_txt
    gltf_texture_info* emission_txt_info = nullptr;

    // reflectance --------------------------
    /// metallic roughnesss
    gltf_material_metallic_roughness* metallic_roughness = nullptr;
    /// specular glossiness
    gltf_material_specular_glossiness* specular_glossiness = nullptr;

    // other textures -----------------------
    /// occlusion texture
    gltf_texture* occlusion_txt = nullptr;
    /// normal texture
    gltf_texture* normal_txt = nullptr;
    /// texture information for collusion_txt
    gltf_texture_info* occlusion_txt_info = nullptr;
    /// texture information for normal_txt
    gltf_texture_info* normal_txt_info = nullptr;

    // other Parameters ---------------------
    /// double sided
    bool double_sided = true;

    // cleanup
    ~gltf_material() {
        if (metallic_roughness) delete metallic_roughness;
        if (specular_glossiness) delete specular_glossiness;
        if (emission_txt_info) delete emission_txt_info;
        if (occlusion_txt_info) delete occlusion_txt_info;
        if (normal_txt_info) delete normal_txt_info;
    }
};

/// Morph information for shapes
struct gltf_shape_morph {
    /// morph position
    std::vector<vec3f> pos;
    /// morph normal
    std::vector<vec3f> norm;
    /// morph tangent
    std::vector<vec3f> tangsp;
    /// default weight (the same for each shape in a mesh)
    float weight = 0;
};

/// Primitives

struct gltf_shape {
    /// name of the mesh that enclosed it
    std::string name = "";
    /// material reference
    gltf_material* mat = nullptr;

    /// vertex position
    std::vector<vec3f> pos;
    /// vertex normal
    std::vector<vec3f> norm;
    /// vertex texcoord
    std::vector<vec2f> texcoord;
    /// vertex additional texcoord
    std::vector<vec2f> texcoord1;
    /// vertex color
    std::vector<vec4f> color;
    /// vertex radius
    std::vector<float> radius;
    /// vertex tangent space
    std::vector<vec4f> tangsp;

    /// vertex skinning weights
    std::vector<vec4f> skin_weights;
    /// vertex skinning joint indices
    std::vector<vec4i> skin_joints;

    /// point elements
    std::vector<int> points;
    /// line elements
    std::vector<vec2i> lines;
    /// triangle elements
    std::vector<vec3i> triangles;

    /// morph targets
    std::vector<gltf_shape_morph*> morph_targets;

    // Cleanup
    ~gltf_shape() {
        for (auto v : morph_targets) delete v;
    }
};

/// Gltf mesh.
struct gltf_mesh {
    /// name
    std::string name = "";
    /// path (only used when writing files on disk with glTF)
    std::string path = "";
    /// primitives references
    std::vector<gltf_shape*> shapes;

    // Cleanup
    ~gltf_mesh() {
        for (auto v : shapes) delete v;
    }
};

// forward declaration
struct gltf_skin;

/// Node in the hierarchy.
struct gltf_node {
    /// name
    std::string name = "";
    /// camera reference
    gltf_camera* cam = nullptr;
    /// mesh reference
    gltf_mesh* msh = nullptr;
    /// mesh reference
    gltf_skin* skn = nullptr;
    /// children
    std::vector<gltf_node*> children;

    /// A floating-point 4x4 transformation matrix stored in column-major order.
    mat4f matrix = identity_mat4f;
    /// The node's unit quaternion rotation in the order (x, y, z, w), where w
    /// is the scalar.
    quat4f rotation = {0, 0, 0, 1};
    /// The node's non-uniform scale.
    vec3f scale = {1, 1, 1};
    /// The node's translation.
    vec3f translation = {0, 0, 0};
    /// morph target weights
    std::vector<float> morph_weights;

    // computed properties ---------------
    /// parent node (computed during update_node_hierarchy())
    gltf_node* parent = nullptr;

    // computed properties ---------------
    /// transform (computed during update_transforms())
    mat4f xform() const { return _xform; }
    /// local transform (computed during update_transforms())
    mat4f local_xform() const { return _local_xform; }
    /// skin transform (computed during update_transforms())
    mat4f skin_xform() const { return _skin_xform; }

    // cached properties -----------------
    // do not access directly: this is likely to change
    /// transform (computed during update_transforms())
    mat4f _xform = identity_mat4f;
    /// local transform (computed during update_transforms())
    mat4f _local_xform = identity_mat4f;
    /// skin transform (computed during update_transforms())
    mat4f _skin_xform = identity_mat4f;

    // cleanup
    ~gltf_node() {
        for (auto v : children) delete v;
    }
};

/// Animation Interpolation
enum struct gltf_animation_interpolation {
    /// linear
    linear = 0,
    /// step function
    step = 1,
    /// catmull-rom spline
    catmull_rom = 2,
    /// cubic bezier spline
    cubic = 3,
};

/// Keyframe data.
struct gltf_animation {
    /// Interpolation
    gltf_animation_interpolation interp = gltf_animation_interpolation::step;
    /// Target nodes
    std::vector<gltf_node*> nodes;
    /// Times
    std::vector<float> time;
    /// Translation
    std::vector<vec3f> translation;
    /// Rotation
    std::vector<quat4f> rotation;
    /// Scale
    std::vector<vec3f> scale;
    /// Weights for morphing
    std::vector<std::vector<float>> morph_weights;
};

/// Animation
struct gltf_animation_group {
    /// Name
    std::string name;
    /// path (only used when writing files on disk with glTF)
    std::string path = "";
    /// Times
    std::vector<gltf_animation*> animations;

    // cleanup
    ~gltf_animation_group() {
        for (auto v : animations) delete v;
    }
};

/// Skin
struct gltf_skin {
    /// name
    std::string name = "";
    /// path (only used when writing files on disk with glTF)
    std::string path = "";
    /// inverse bind matrix
    std::vector<mat4f> pose_matrices;
    /// joints
    std::vector<gltf_node*> joints;
    /// skeleton root node
    gltf_node* root = nullptr;
};

/// Gltf scene
struct gltf_scene {
    /// name
    std::string name = "";
    /// instances
    std::vector<gltf_node*> nodes;
};

/// Gltf model. Objects are shared between scenes.
/// Scenes and nodes are missing for mesh-only assets.
struct gltf_scene_group {
    /// default scene (null if not present)
    gltf_scene* default_scene = nullptr;
    /// cameras
    std::vector<gltf_camera*> cameras;
    /// materials
    std::vector<gltf_material*> materials;
    /// textures
    std::vector<gltf_texture*> textures;
    /// meshes
    std::vector<gltf_mesh*> meshes;
    /// scenes
    std::vector<gltf_scene*> scenes;
    /// nodes
    std::vector<gltf_node*> nodes;
    /// nodes
    std::vector<gltf_animation_group*> animations;
    /// skins
    std::vector<gltf_skin*> skins;

    // cleanup
    ~gltf_scene_group() {
        for (auto v : cameras) delete v;
        for (auto v : materials) delete v;
        for (auto v : textures) delete v;
        for (auto v : meshes) delete v;
        for (auto v : scenes) delete v;
        for (auto v : nodes) delete v;
        for (auto v : animations) delete v;
        for (auto v : skins) delete v;
    }
};

/// Load scene
///
/// - Parameters:
///     - filename: filename
///     - load_textures: whether to load textures (default to false)
///     - skip_missing: whether to skip missing buffers and textures
/// - Returns:
///     - scene (nullptr on error)
gltf_scene_group* load_scenes(
    const std::string& filename, bool load_textures, bool skip_missing = true);

/// Save scene
///
/// - Parameters:
///     - filename: filename
///     - buffer_uri: name of the main buffer
///     - scn: scene data to save
///     - save_textures: whether to save textures (default to false)
///     - separate_buffers: save separate buffers for each mesh
void save_scenes(const std::string& filename, const std::string& buffer_uri,
    const gltf_scene_group* scn, bool save_textures,
    bool separate_buffers = false);

/// Update node hierarchy
void update_node_hierarchy(gltf_scene_group* scn);

/// Update node trasforms
void update_transforms(gltf_scene_group* scn);

/// Update animated node
void update_animated_transforms(gltf_scene_group* scns, float time);

/// Get a list of nodes with meshes
std::vector<gltf_node*> get_mesh_nodes(const gltf_scene* scn);

/// Get a list of nodes with cameras
std::vector<gltf_node*> get_camera_nodes(const gltf_scene* scn);

/// Animation times
vec2f get_animation_bounds(const gltf_scene_group* scn);

/// Skin transforms (local-to-object) from the node transform that instances the
/// skin
std::vector<mat4f> get_skin_transforms(const gltf_skin* sk, const mat4f& xform);

/// Compute shape morphing
void compute_morphing_deformation(const gltf_shape* shp,
    const std::vector<float>& weights, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec4f>& tangsp);

/// Computes a scene bounding box
bbox3f compute_scene_bounds(const gltf_scene_group* scn);

/// Add missing data to the scene.
void add_normals(gltf_scene_group* scn);

/// Add missing data to the scene.
void add_radius(gltf_scene_group* scn, float radius);

/// Add missing data to the scene.
void add_tangent_space(gltf_scene_group* scn);

/// Add missing data to the scene.
void add_nodes(gltf_scene_group* scn);

/// Add missing data to the scene.
void add_scene(gltf_scene_group* scn);

/// Add missing data to the scene.
void add_texture_data(gltf_scene_group* scn);

/// Add missing data to the scene.
void add_names(gltf_scene_group* scn);

/// Add a default camera that views the entire scene.
void add_default_cameras(gltf_scene_group* scn);

/// Set unique path names for outputting separate buffers
void add_unique_path_names(
    gltf_scene_group* scns, const std::string& buffer_uri);

/// Convert materials to spec gloss
void add_spec_gloss(gltf_scene_group* scns);

/// Convert a gltf asset to flattened group of scene.
gltf_scene_group* gltf_to_scenes(const glTF* gltf, int scene_idx = -1);

/// Convert a flattened group of scene into a gltf. If separate_buffers,
/// creates a separate buffer for each each and animation and
/// prepend buffer_uri to its name.
glTF* scenes_to_gltf(const gltf_scene_group* fl_gltf,
    const std::string& buffer_uri, bool separate_buffers = false);

/// Validate a gltf. Missing many validation as of this version.
std::vector<std::pair<std::string, std::string>> validate_gltf(glTF* gltf);

}  // namespace ygl

#endif
