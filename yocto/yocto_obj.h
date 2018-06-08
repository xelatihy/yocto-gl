//
// # Yocto/Obj: Tiny C++ Library for OBJ loading/saving
//
// Yocto/Obj is a library for loading and saving OBJ files for 3D applications,
// support for points, lines, triangles and general polygons and all materials
// properties. Yocto/Obj defined also a few extensions to easily create demos
// such as per-vertex color and radius, cameras, environment maps and instances.
// Can use either a low-level OBJ representation, from this files,
// or a high level flattened representation included in Yocto/Scene.
//
// Both in reading and writing, OBJ has no clear convention on the orientation
// of textures Y axis. So in many cases textures appears flipped. To handle
// that, use the option to flip textures coordinates on either saving or
// loading. By default texture coordinates are flipped since this seems
// the convention found on test cases collected on the web. The value Tr
// has similar problems, since its relation to opacity is software specific.
// Again we let the user choose the conversion and set the default to the
// one found on the web.
//
// Yocto/Obj stores all parsed data in its data structures. This design
// decision does not scale particularly well with very large models. The main
// reasons is the OBJ file format vertices are global for the entire scene and
// not single objects. To avoid memory allocation problems, we internally store
// vertices in queues during reading and make them unique per object right after
// parsing.
//
// In the high level interface, shapes are indexed meshes and are described
// by arrays of vertex indices for points/lines/triangles and arrays for vertex
// positions, normals, texcoords, color and radius. The latter two as
// extensions. Since OBJ is a complex formats that does not match well with
// current GPU rendering / path tracing algorithms, we adopt a simplification
// similar to other single file libraries:
// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
//   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
//   vertex duplication happens thought for same triplets
// 2. we split shapes on changes to groups and materials, instead of keeping
//   per-face group/material data; this makes the data usable right away in
//   a GPU viewer; this is not a major limitation if we accept the previous
//   point that already changes shapes topology.
//
// ## Usage
//
// 1. load a obj data with `load_obj()`; can load also textues
// 2. look at the `obj_XXX` data structures for access to individual elements
// 3. use obj back to disk with `save_obj()`; can also save textures
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

#ifndef _YGL_OBJ_H_
#define _YGL_OBJ_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

// enable fast OBJ parsing
#ifndef YGL_FASTOBJ
#define YGL_FASTOBJ 1
#endif

#include "yocto_math.h"

#include <deque>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// WAVEFRONT OBJ SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

// Obj element type.
enum struct obj_element_type : uint8_t {
    point = 1,
    line = 2,
    face = 3,
    bezier = 4
};

// Obj element (point/polyline/polygon)
struct obj_element {
    uint32_t start;          // starting vertex index
    uint8_t size;            // number of vertices
    obj_element_type type;   // element type
    uint16_t material = 0;   // material id
    uint16_t group = 0;      // group id
    uint16_t smoothing = 0;  // smoothing group id
};

// Obj group properties.
struct obj_group_props {
    std::string name = "";     // group name
    std::string matname = "";  // material name
    bool faceted = false;      // faceted or smooth
};

// Obj object.
struct obj_object {
    std::string name;                    // name
    std::vector<std::string> materials;  // materials
    std::vector<std::string> groups;     // groups
    std::vector<bool> smoothing;         // smoothing groups
    std::vector<int> verts_pos;          // element pos
    std::vector<int> verts_norm;         // element norm
    std::vector<int> verts_texcoord;     // element texcoord
    std::vector<obj_element> elems;      // elements
    frame3f frame = identity_frame3f;    // frame [extension]
    vec2i subdiv = zero2i;  // type/level of subdivision [extension]
    // Properties not explicitly handled [extension].
    std::unordered_map<std::string, std::vector<float>> props;
};

// Obj texture. Texture data is loaded only if desired.
struct obj_texture {
    std::string path;        // path
    int width = 0;           // width
    int height = 0;          // height
    std::vector<vec4f> img;  // image data
};

// Obj texture information.
struct obj_texture_info {
    std::string path = "";  // file path
    bool clamp = false;     // clamp to edge
    float scale = 1;        // scale for bump/displacement
    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<float>> props;
};

// Comparison for texture info.
inline bool operator==(const obj_texture_info& a, const obj_texture_info& b) {
    if (a.path.empty() && b.path.empty()) return true;
    if (a.path != b.path) return false;
    return a.clamp == b.clamp && a.scale == b.scale && a.props == b.props;
}

// Obj material.
struct obj_material {
    std::string name;  // name
    int illum = 0;     // MTL illum mode

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f ka = {0, 0, 0};  // ambient color
    vec3f kd = {0, 0, 0};  // diffuse color
    vec3f ks = {0, 0, 0};  // specular color
    vec3f kr = {0, 0, 0};  // reflection color
    vec3f kt = {0, 0, 0};  // transmission color
    float ns = 0;          // Phong exponent color
    float ior = 1;         // index of refraction
    float op = 1;          // opacity
    float rs = -1;         // roughness (-1 not defined)
    float km = -1;         // metallic  (-1 not defined)

    obj_texture_info ke_txt;    // emission texture
    obj_texture_info ka_txt;    // ambient texture
    obj_texture_info kd_txt;    // diffuse texture
    obj_texture_info ks_txt;    // specular texture
    obj_texture_info kr_txt;    // reflection texture
    obj_texture_info kt_txt;    // transmission texture
    obj_texture_info ns_txt;    // Phong exponent texture
    obj_texture_info op_txt;    // opacity texture
    obj_texture_info rs_txt;    // roughness texture
    obj_texture_info km_txt;    // metallic texture
    obj_texture_info ior_txt;   // ior texture
    obj_texture_info occ_txt;   // occlusion map
    obj_texture_info bump_txt;  // bump map
    obj_texture_info disp_txt;  // displacement map
    obj_texture_info norm_txt;  // normal map

    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Obj camera [extension].
struct obj_camera {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    bool ortho = false;                // orthographic
    float width = 0.036f;              // film width (default to 35mm)
    float height = 0.024f;             // film height (default to 35mm)
    float focal = 0.050f;              // focal length
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float aperture = 0;                // lens aperture
    float focus = flt_max;             // focus distance
};

// Obj environment [extension].
struct obj_environment {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    vec3f ke = zero3f;                 // emission color
    obj_texture_info ke_txt;           // emission texture
};

// Obj scene.
struct obj_scene {
    std::deque<vec3f> pos;       // vertex positions
    std::deque<vec3f> norm;      // vertex normals
    std::deque<vec2f> texcoord;  // vertex texcoords

    std::vector<obj_object*> objects;            // objects
    std::vector<obj_material*> materials;        // materials
    std::vector<obj_texture*> textures;          // textures
    std::vector<obj_camera*> cameras;            // cameras [extension]
    std::vector<obj_environment*> environments;  // environments [extension]

    // Cleanup.
    ~obj_scene() {
        for (auto v : objects) delete v;
        for (auto v : materials) delete v;
        for (auto v : textures) delete v;
        for (auto v : cameras) delete v;
        for (auto v : environments) delete v;
    }
};

// Load an OBJ from file `filename`. Split shapes at material and group
// boundaries, if `split_shapes` is true.
// Load textures if `load_textures` is true, and report errors only if
// `skip_missing` is false. Texture coordinates and material Tr are flipped
// if `flip_texcoord` and `flip_tp` are respectively true.
obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool flip_texcoord = true, bool flip_tr = true);

// Save an OBJ to file `filename`. Save textures if `save_textures` is true,
// and report errors only if `skip_missing` is false.
// Texture coordinates and material Tr are flipped if `flip_texcoord` and
// `flip_tp` are respectively true.
void save_obj(const std::string& filename, const obj_scene* obj,
    bool flip_texcoord = true, bool flip_tr = true);

// Load OBJ texture images.
void load_obj_textures(
    obj_scene* obj, const std::string& dirname, bool skip_missing = true);
// Save OBJ texture images.
void save_obj_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing = true);

// Callback object
struct obj_callbacks {
    bool (*object)(void*, const std::string&) = nullptr;  // object callback
    bool (*group)(void*, const std::string&) = nullptr;   // group callback
    bool (*usemat)(void*, const std::string&) = nullptr;  // use material cb
    bool (*smoothing)(void*, bool) = nullptr;             // smoothing callback
    bool (*vertex)(void*, vec3f) = nullptr;               // vertex callback
    bool (*normal)(void*, vec3f) = nullptr;               // normal callback
    bool (*texcoord)(void*, vec2f) = nullptr;             // texcoord callback
    bool (*face)(void*, int, const vec3i*) = nullptr;     // face callback
    bool (*line)(void*, int, const vec3i*) = nullptr;     // line callback
    bool (*point)(void*, int, const vec3i*) = nullptr;    // point callback
    bool (*material)(
        void*, const obj_material&) = nullptr;           // material callback
    bool (*camera)(void*, const obj_camera&) = nullptr;  // camera callback
    bool (*environment)(void*, const obj_environment&) = nullptr;  // envmap cb
};

// Load OBJ with callbacks
void load_obj(
    const std::string& filename, const obj_callbacks& callbacks, void* ctx);

}  // namespace ygl

#endif
