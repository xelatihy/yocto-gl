///
/// # Yocto/Obj
///
/// Wavefront OBJ/MTL loader and writer with support for points,
/// lines, triangles and general polygons and all materials properties.
/// Contains also a few extensions to easily create demos such as per-vertex
/// color and radius, cameras, environment maps and instances.
/// Can use either a low-level OBJ representation or a high level flattened
/// representation.
///
/// Both in reading and writing, OBJ has no clear convention on the orientation
/// of textures Y axis. So in many cases textures appears flipped. To handle
/// that, use the option to flip textures coordinates on either saving or
/// loading. By default texture coordinates are flipped since this seems
/// the convention found on test cases collected on the web. The value Tr
/// has similar problems, since its relation to opacity is software specific.
/// Again we let the user chose the convension and set the default to the
/// one found on the web.
///
/// In the high level interface, shapes are indexed meshes and are described
/// by arrays of vertex indices for points/lines/triangles and arrays for vertex
/// positions, normals, texcoords, color and radius. The latter two as
/// extensions. Since OBJ is a complex formats that does not match well with
/// current GPU rendering / path tracing algorithms, we adopt a simplification
/// similar to other single file libraries:
/// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
///   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
///   vertex dusplication happens thought for same triplets
/// 2. we split shapes on changes to groups and materials, instead of keeping
///   per-face group/material data; this makes the data usable right away in
///   a GPU viewer; this is not a major limitation if we accept the previous
///   point that already changes shapes topology.
///
/// This library depends in yocto_math.h. Texture loading depends on
/// yocto_image. If the texture loading dependency is not desired, it can be
/// disabled by defining YOBJ_NO_IMAGE before including this file.
///
///
/// ## Usage Of High-Level Interface
///
/// 1. load a scene with `load_scene()`
/// 2. look at the `scene` data structures for access to individual elements
/// 3. can also manipulate the scene by adding missing data with `add_XXX()`
///    functions
/// 4. since OBJ does natively support mesh transfotms, which we support with
///    instances, use `flatten_instaces()` or `add_instances()` to go back
///    and fourth
/// 5. use `save_scene()` ti write the data to disk
///
/// ## Usage Of Low-Level Interface
///
/// 1. load a obj data with `load_obj()`; can load also textues
/// 2. look at the `obj_XXX` data structures for access to individual elements
/// 3. use obj back to disk with `save_obj()`; can also save textures
/// 4. conversion from low- to -high-level data structures with
///    `scene_to_obj()` and `obj_to_scene()`
///
///
/// ## History
///
/// - v 0.30: support for smoothing groups
/// - v 0.29: use reference interface for textures
/// - v 0.28: add function to split meshes into single shapes
/// - v 0.27: explicit transforms
/// - v 0.26: added interpreting of illum in scene conversions
/// - v 0.25: added convention for Tr
/// - v 0.24: remove exception from code and add explicit error handling
/// - v 0.23: texture have always 4 channels
/// - v 0.22: change variable names for compilation on gcc
/// - v 0.21: bug fixes
/// - v 0.20: use yocto_math in the interface and remove inline compilation
/// - v 0.19: add missing bounding box computation and missing data functions
/// - v 0.18: prioritize high-level interface
/// - v 0.17: name cleanup in both interface to better align with glTF
/// - v 0.16: change flattened data structure to use pointers
/// - v 0.15: added unknown properties std::map for materials
/// - v 0.14: added extension to store tetrahedral meshes
/// - v 0.13: started adding physics extension to materials
/// - v 0.12: change texture loading by flipping uvs rather than images
/// - v 0.11: use yocto_image for texture handling.
/// - v 0.10: switch to .h/.cpp pair
/// - v 0.9: bug fixes and optionally texture skipping
/// - v 0.8: high level interface uses grouping
/// - v 0.7: doxygen comments
/// - v 0.6: bug fixes
/// - v 0.5: removed options to force image formats (image library not reliable)
/// - v 0.4: [major API change] move to modern C++ interface
/// - v 0.3: new API internals and C++ interface
/// - v 0.2: removal of C interface
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace yobj {}

//
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
//

#ifndef _YOBJ_H_
#define _YOBJ_H_

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// Reading and Writing support for Wavefront OBJ.
///
namespace yobj {

// -----------------------------------------------------------------------------
// HIGH-LEVEL INTERFACE
// -----------------------------------------------------------------------------

///
/// Property map
///
template <typename T>
using property_map = std::map<std::string, std::vector<T>>;

///
/// Scene Texture
///
struct texture {
    /// path
    std::string path;
    /// if loaded, ldr image
    ym::image4b ldr;
    /// if loaded, hdr image
    ym::image4f hdr;

    /// get texture width
    int width() const {
        if (ldr) return ldr.width();
        if (hdr) return hdr.width();
        return 0;
    }
    /// get texture height
    int height() const {
        if (ldr) return ldr.height();
        if (hdr) return hdr.height();
        return 0;
    }
};

///
/// Scene Texture Additional Information
///
struct texture_info {
    /// clamping
    bool clamp = false;
    /// bump scale
    float bump_scale = 1;
    /// unknown string props
    property_map<std::string> unknown_props;
};

///
/// Scene Material
///
struct material {
    // whole material data -------------------
    /// material name
    std::string name;

    // color information ---------------------
    /// emission color
    ym::vec3f ke = {0, 0, 0};
    /// diffuse color
    ym::vec3f kd = {0, 0, 0};
    /// specular color
    ym::vec3f ks = {0, 0, 0};
    /// transmission color
    ym::vec3f kt = {0, 0, 0};
    /// roughness
    float rs = 0.0001;
    /// opacity
    float opacity = 1;

    // textures -------------------------------
    /// emission texture
    texture* ke_txt = nullptr;
    /// diffuse texture
    texture* kd_txt = nullptr;
    /// specular texture
    texture* ks_txt = nullptr;
    /// transmission texture
    texture* kt_txt = nullptr;
    /// roughness texture
    texture* rs_txt = nullptr;
    /// opacity texture
    texture* op_txt = nullptr;
    /// bump map texture (heighfield)
    texture* bump_txt = nullptr;
    /// displacement map texture (heighfield)
    texture* disp_txt = nullptr;
    /// normal texture
    texture* norm_txt = nullptr;

    // texture information ---------------------
    /// emission texture
    texture_info ke_txt_info = {};
    /// diffuse texture
    texture_info kd_txt_info = {};
    /// specular texture
    texture_info ks_txt_info = {};
    /// transmission texture
    texture_info kt_txt_info = {};
    /// roughness texture
    texture_info rs_txt_info = {};
    /// bump map texture (heighfield)
    texture_info bump_txt_info = {};
    /// displacement map texture (heighfield)
    texture_info disp_txt_info = {};
    /// normal texture
    texture_info norm_txt_info = {};

    // unknown properties ---------------------
    /// unknown string props
    property_map<std::string> unknown_props;
};

///
/// Shape. May contain only one of the points/lines/triangles.
///
struct shape {
    /// name of the group that enclosed it
    std::string name = "";
    /// material
    material* mat = nullptr;

    // shape elements -------------------------
    /// points
    std::vector<int> points;
    /// lines
    std::vector<ym::vec2i> lines;
    /// triangles
    std::vector<ym::vec3i> triangles;
    /// tetrahedrons
    std::vector<ym::vec4i> tetras;

    // vertex data ----------------------------
    /// per-vertex position (3 float)
    std::vector<ym::vec3f> pos;
    /// per-vertex normals (3 float)
    std::vector<ym::vec3f> norm;
    /// per-vertex texcoord (2 float)
    std::vector<ym::vec2f> texcoord;
    /// [extension] per-vertex color (4 float)
    std::vector<ym::vec4f> color;
    /// [extension] per-vertex radius (1 float)
    std::vector<float> radius;
    /// [extension] per-vertex tangent space (4 float)
    std::vector<ym::vec4f> tangsp;
};

///
/// Mesh
///
struct mesh {
    // name
    std::string name;
    /// primitives
    std::vector<shape*> shapes;

    /// cleanup
    ~mesh();
};

///
/// Mesh instance.
///
struct instance {
    // name
    std::string name;
    /// translation
    ym::vec3f translation = {0, 0, 0};
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// scale
    ym::vec3f scale = {1, 1, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;
    /// mesh instances
    mesh* msh = nullptr;

    /// xform
    ym::mat4f xform() const {
        return ym::translation_mat4(translation) * ym::rotation_mat4(rotation) *
               ym::scaling_mat4(scale) * matrix;
    }
};

///
/// Scene Camera
///
struct camera {
    /// name
    std::string name;
    /// translation
    ym::vec3f translation = {0, 0, 0};
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;
    /// ortho cam
    bool ortho = false;
    /// vertical field of view
    float yfov = 2;
    /// aspect ratio
    float aspect = 16.0f / 9.0f;
    /// focus distance
    float focus = 1;
    /// lens aperture
    float aperture = 0;

    /// xform
    ym::mat4f xform() const {
        return ym::translation_mat4(translation) * ym::rotation_mat4(rotation) *
               matrix;
    }
};

///
/// Envinonment map
///
struct environment {
    /// name
    std::string name;
    /// index of material in material array
    material* mat = nullptr;
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;

    /// xform
    ym::mat4f xform() const { return ym::rotation_mat4(rotation) * matrix; }
};

///
/// Scene
///
struct scene {
    /// shape array
    std::vector<mesh*> meshes;
    /// instance array
    std::vector<instance*> instances;
    /// material array
    std::vector<material*> materials;
    /// texture array
    std::vector<texture*> textures;
    /// camera array
    std::vector<camera*> cameras;
    /// environment array
    std::vector<environment*> environments;

    /// cleanup
    ~scene();
};

///
/// Load scene
///
/// - Parameters:
///     - filename: filename
///     - load_textures: whether to load textures (default to false)
///     - flip_texcoord: whether to flip the v coordinate
///     - facet_non_smooth: duplicate vertices if smoothing off
///     - skip_missing: skip missing textures
///     - flip_tr: whether to flip tr
///     - err: if set, store error message on error
/// - Returns:
///     - scene (nullptr on error)
///
scene* load_scene(const std::string& filename, bool load_textures,
    bool skip_missing = true, bool flip_texcoord = true,
    bool facent_non_smooth = false, bool flip_tr = true,
    std::string* err = nullptr);

///
/// Save scene
///
/// - Parameters:
///     - filename: filename
///     - scn: scene data to save
///     - save_textures: whether to save textures (default to false)
///     - flip_texcoord: whether to flip the v coordinate
///     - flip_tr: whether to flip tr
///     - err: if set, store error message on error
/// - Returns:
///     - whether an error occurred
///
bool save_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool flip_texcoord = true, bool flip_tr = true,
    std::string* err = nullptr);

#ifndef YOBJ_NO_IMAGE

///
/// Loads textures for an scene.
///
/// - Parameters:
///     - scn: scene to load textures into
///     - dirname: base directory name for texture files
///     - skip_missing: whether to skip missing textures or stops with error
///     - err: if set, store error message on error
/// - Returns:
///     - whether an error occurred
///
bool load_textures(scene* scn, const std::string& dirname,
    bool skip_missing = false, std::string* err = nullptr);

///
/// Saves textures for an scene.
///
/// - Parameters:
///     - scn: scene to write textures from
///     - dirname: base directory name for texture files
///     - skip_missing: whether to skip missing textures or stops with error
///     - err: if set, store error message on error
/// - Returns:
///     - whether an error occurred
///
bool save_textures(const scene* scn, const std::string& dirname,
    bool skip_missing = false, std::string* err = nullptr);

#endif

///
/// Computes a scene bounding box
///
ym::bbox3f compute_scene_bounds(const scene* scn);

///
/// Add missing data to the scene.
///
void add_normals(scene* scn);

///
/// Add missing data to the scene.
///
void add_radius(scene* scn, float radius);

///
/// Add missing data to the scene.
///
void add_tangent_space(scene* scn);

///
/// Add missing data to the scene.
///
void add_texture_data(scene* scn);

///
/// Add missing data to the scene.
///
void add_instances(scene* scn);

///
/// Add missing data to the scene.
///
void add_names(scene* scn);

///
/// Add a default camera that views the entire scene.
///
void add_default_camera(scene* scn);

///
/// Flatten scene instances into separate meshes.
///
void flatten_instances(scene* scn);

///
/// Split meshes into single shapes
///
void split_shapes(scene* scn);

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------

///
/// Face vertex
///
struct obj_vertex {
    /// position
    int pos;
    /// texcoord
    int texcoord;
    /// normal
    int norm;
    /// color [extension]
    int color;
    /// radius [extension]
    int radius;

    /// Constructor (copies members initializing missing ones to -1)
    obj_vertex(int pos = -1, int texcoord = -1, int norm = -1, int color = -1,
        int radius = -1)
        : pos(pos)
        , texcoord(texcoord)
        , norm(norm)
        , color(color)
        , radius(radius) {}
};

///
/// element type
///
enum struct obj_element_type : uint16_t {
    /// lists of points
    point = 1,
    /// polylines
    line = 2,
    /// polygon faces
    face = 3,
    /// tetrahedrons
    tetra = 4,
};

///
/// Element vertex indices
///
struct obj_element {
    /// starting vertex index
    uint32_t start;
    /// element type
    obj_element_type type;
    /// number of vertices
    uint16_t size;
};

///
/// Element group
///
struct obj_group {
    // group data ---------------------------
    /// material name
    std::string matname;
    /// group name
    std::string groupname;
    /// smoothing
    bool smoothing = true;

    // element data -------------------------
    /// element vertices
    std::vector<obj_vertex> verts;
    /// element faces
    std::vector<obj_element> elems;
};

///
/// Obj object
///
struct obj_object {
    // object data --------------------------
    /// object name
    std::string name;

    // element data -------------------------
    /// element groups
    std::vector<obj_group> groups;
};

///
/// OBJ material
///
struct obj_material {
    // whole material data ------------------
    /// material name
    std::string name;
    /// MTL illum mode
    int illum = 0;

    // color information --------------------
    /// emission color
    ym::vec3f ke = {0, 0, 0};
    /// ambient color
    ym::vec3f ka = {0, 0, 0};
    /// diffuse color
    ym::vec3f kd = {0, 0, 0};
    /// specular color
    ym::vec3f ks = {0, 0, 0};
    /// reflection color
    ym::vec3f kr = {0, 0, 0};
    /// transmision color
    ym::vec3f kt = {0, 0, 0};
    /// phong exponent for ks
    float ns = 1;
    /// index of refraction
    float ior = 1;
    /// opacity
    float op = 1;

    // texture names for the above properties
    /// emission texture
    std::string ke_txt;
    /// ambient texture
    std::string ka_txt;
    /// diffuse texture
    std::string kd_txt;
    /// specular texture
    std::string ks_txt;
    /// reflection texture
    std::string kr_txt;
    /// transmission texture
    std::string kt_txt;
    /// specular exponent texture
    std::string ns_txt;
    /// opacity texture
    std::string op_txt;
    /// index of refraction
    std::string ior_txt;
    /// bump map texture (heighfield)
    std::string bump_txt;
    /// displacement map texture (heighfield)
    std::string disp_txt;
    /// normal map texture
    std::string norm_txt;

    // texture information ---------------------
    /// emission texture
    property_map<std::string> ke_txt_info = {};
    /// ambient texture
    property_map<std::string> ka_txt_info = {};
    /// diffuse texture
    property_map<std::string> kd_txt_info = {};
    /// specular texture
    property_map<std::string> ks_txt_info = {};
    /// reflection texture
    property_map<std::string> kr_txt_info = {};
    /// transmission texture
    property_map<std::string> kt_txt_info = {};
    /// specular exponent texture
    property_map<std::string> ns_txt_info = {};
    /// opacity texture
    property_map<std::string> op_txt_info = {};
    /// index of refraction
    property_map<std::string> ior_txt_info = {};
    /// bump map texture (heighfield)
    property_map<std::string> bump_txt_info = {};
    /// displacement map texture (heighfield)
    property_map<std::string> disp_txt_info = {};
    /// normal texture
    property_map<std::string> norm_txt_info = {};

    // unknown properties ---------------------
    /// unknown string props
    property_map<std::string> unknown_props;
};

///
/// Camera [extension]
///
struct obj_camera {
    /// camera name
    std::string name;
    /// translation
    ym::vec3f translation = {0, 0, 0};
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// scale
    ym::vec3f scale = {1, 1, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;
    /// orthografic camera
    bool ortho = false;
    /// vertical field of view
    float yfov = 2 * std::atan(0.5f);
    /// aspect ratio
    float aspect = 16.0f / 9.0f;
    /// lens aperture
    float aperture = 0;
    /// focus distance
    float focus = 1;
};

///
/// Environment [extension]
///
struct obj_environment {
    /// environment name
    std::string name;
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;
    /// material name
    std::string matname;
};

///
/// Instance [extension]
///
struct obj_instance {
    /// instance name
    std::string name;
    /// translation
    ym::vec3f translation = {0, 0, 0};
    /// rotation
    ym::quat4f rotation = {0, 0, 0, 1};
    /// scale
    ym::vec3f scale = {1, 1, 1};
    /// generic transform matrix
    ym::mat4f matrix = ym::identity_mat4f;
    /// object name
    std::string meshname;
};

///
/// OBJ asset
///
struct obj {
    // vertex data -------------------------
    /// vertex positions
    std::vector<ym::vec3f> pos;
    /// vertex normals
    std::vector<ym::vec3f> norm;
    /// vertex texcoord
    std::vector<ym::vec2f> texcoord;
    /// vertex color [extension]
    std::vector<ym::vec4f> color;
    /// vertex radius [extension]
    std::vector<float> radius;

    // scene objects -----------------------
    /// objects
    std::vector<obj_object> objects;
    /// materials
    std::vector<obj_material> materials;
    /// cameras [extension]
    std::vector<obj_camera> cameras;
    /// env maps [extension]
    std::vector<obj_environment> environments;
    /// instances [extension]
    std::vector<obj_instance> instances;
};

///
/// Load OBJ
///
/// - Parameters:
///     - filename: filename
///     - flip_texcoord: whether to flip the v coordinate
///     - flip_tr: whether to flip the Tr value
///     - err: if set, store error message on error
/// - Return:
///     - obj (nullptr on error)
///
obj* load_obj(const std::string& filename, bool flip_texcoord = true,
    bool flip_tr = true, std::string* err = nullptr);

///
/// Load MTL
///
/// - Parameters:
///     - filename: filename
///     - flip_tr: whether to flip the Tr value
///     - err: if set, store error message on error
/// - Return:
///     - loaded materials (empty on error)
///
std::vector<obj_material> load_mtl(const std::string& filename,
    bool flip_tr = true, std::string* err = nullptr);

///
/// Save OBJ
///
/// - Parameters:
///     - filename: filename
///     - asset: obj data to save
///     - flip_texcoord: whether to flip the v coordinate
///     - flip_tr: whether to flip the Tr value
///     - err: if set, store error message on error
/// - Returns:
///     - whether an error occurred
///
bool save_obj(const std::string& filename, const obj* asset,
    bool flip_texcoord = true, bool flip_tr = true, std::string* err = nullptr);

///
/// Save MTL (@deprecated interface)
///
bool save_mtl(const std::string& filename,
    const std::vector<obj_material>& materials, bool flip_tr = true,
    std::string* err = nullptr);

///
/// Converts an OBJ into a scene.
///
/// - Parameters:
///     - obj: obj to be flattened
///     - facet_non_smooth: duplicate vertices if smoothing off
///
scene* obj_to_scene(const obj* obj, bool facet_non_smooth);

///
/// Save a scene in an OBJ file.
///
/// - Parameters:
///     - scn: scene to convert
///
obj* scene_to_obj(const scene* scn);

}  // namespace yobj

#endif
