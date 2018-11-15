//
// # Yocto/SceneIO: Tiny library for Yocto/Scene input and output
//
// Yocto/SceneIO provides loading and saving functionality for scenes
// in Yocto/GL. We support a simple to use JSON format, PLY, OBJ and glTF.
// The JSON serialization is a straight copy of the in-memory scene data.
// To speed up testing, we also support a binary format that is a dump of
// the current scene. This format should not be use for archival though.
//
// We do not use exception as the API for returning errors, although they might
// be used internally in the implementastion of the methods. In load functions,
// as error is signaled by returning an empty object or a null pointer. In
// save functions, errors are returned with the supplied boolean. In the future,
// we will also provide return types with error codes.
//
// ## Scene Loading and Saving
//
// 1. load a scene with `load_json_scene()` and save it with `save_json_scene()`
// 2. load and save OBJs with `load_obj_scene()` and `save_obj_scene()`
// 3. load and save glTFs with `load_gltf_scene()` and `save_gltf_scene()`
// 4. load and save binary dumps with `load_ybin_scene()`and `save_ybin_scene()`
// 5. if desired, the function `load_scene()` and `save_scene()` will either
//    load using the internal format or convert on the fly using on the
//    supported conversions
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

#ifndef _YOCTO_SCENEIO_H_
#define _YOCTO_SCENEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Scene load options
struct load_scene_options {
    bool          skip_textures             = false;
    bool          exit_on_error             = false;
    bool          obj_split_shapes          = true;
    bool          obj_preserve_face_varying = false;
    bool          assign_texture_opacity    = true;
    atomic<bool>* cancel_flag               = nullptr;
    bool          run_serially              = false;
};
// Scene save options
struct save_scene_options {
    bool          skip_textures = false;
    bool          exit_on_error = false;
    atomic<bool>* cancel_flag   = nullptr;
    bool          run_serially  = false;
};

// Load/save a scene in the supported formats.
bool load_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene in the builtin JSON format.
bool load_json_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_json_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to OBJ.
bool load_obj_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_obj_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to glTF.
bool load_gltf_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_gltf_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
bool load_pbrt_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_pbrt_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a binary dump useful for very fast scene IO. This format is not
// an archival format and should only be used as an intermediate format.
bool load_ybin_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_ybin_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// MESH IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Load/Save a mesh
bool load_mesh(const string& filename, vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool force_triangles);
bool save_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii = false);

// Load/Save a ply mesh
bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<vec4f>& color, vector<float>& radius,
    bool force_triangles);
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, bool force_triangles,
    bool flip_texcoord = true);
bool save_obj_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord = true);

// Load/Save a face-varying mesh
bool load_facevarying_mesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<int>& quads_materials);
bool save_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_facevarying_mesh(const string& filename,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    vector<int>& quads_materials, bool flip_texcoord = true);
bool save_obj_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool flip_texcoord = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// OBJ vertex
struct obj_vertex {
    int position     = 0;
    int texturecoord = 0;
    int normal       = 0;
};

// Obj texture information.
struct obj_texture_info {
    string path  = "";     // file path
    bool   clamp = false;  // clamp to edge
    float  scale = 1;      // scale for bump/displacement
    // Properties not explicitly handled.
    unordered_map<string, vector<float>> props;
};

// Obj material.
struct obj_material {
    string name;       // name
    int    illum = 0;  // MTL illum mode

    // base values
    vec3f ke  = {0, 0, 0};  // emission color
    vec3f ka  = {0, 0, 0};  // ambient color
    vec3f kd  = {0, 0, 0};  // diffuse color
    vec3f ks  = {0, 0, 0};  // specular color
    vec3f kr  = {0, 0, 0};  // reflection color
    vec3f kt  = {0, 0, 0};  // transmission color
    float ns  = 0;          // Phong exponent color
    float ior = 1;          // index of refraction
    float op  = 1;          // opacity
    float rs  = -1;         // roughness (-1 not defined)
    float km  = -1;         // metallic  (-1 not defined)

    // textures
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

    // volume data [extension]
    vec3f ve = {0, 0, 0};  // volume emission
    vec3f va = {0, 0, 0};  // albedo: scattering / (absorption + scattering)
    vec3f vd = {0, 0, 0};  // density: absorption + scattering
    float vg = 0;          // phase function shape

    // volume textures [extension]
    obj_texture_info vd_txt;  // density

    // Properties not explicitly handled.
    unordered_map<string, vector<string>> props;
};

// Obj camera [extension].
struct obj_camera {
    string  name;                         // name
    frame3f frame    = identity_frame3f;  // transform
    bool    ortho    = false;             // orthographic
    float   width    = 0.036f;            // film size (default to 35mm)
    float   height   = 0.024f;            // film size (default to 35mm)
    float   focal    = 0.050f;            // focal length
    float   aspect   = 16.0f / 9.0f;      // aspect ratio
    float   aperture = 0;                 // lens aperture
    float   focus    = float_max;         // focus distance
};

// Obj environment [extension].
struct obj_environment {
    string           name;                      // name
    frame3f          frame = identity_frame3f;  // transform
    vec3f            ke    = zero3f;            // emission color
    obj_texture_info ke_txt;                    // emission texture
};

// Obj callbacks
struct obj_callbacks {
    function<void(const vec3f&)>              vert        = {};
    function<void(const vec3f&)>              norm        = {};
    function<void(const vec2f&)>              texcoord    = {};
    function<void(const vector<obj_vertex>&)> face        = {};
    function<void(const vector<obj_vertex>&)> line        = {};
    function<void(const vector<obj_vertex>&)> point       = {};
    function<void(const string& name)>        object      = {};
    function<void(const string& name)>        group       = {};
    function<void(const string& name)>        usemtl      = {};
    function<void(const string& name)>        smoothing   = {};
    function<void(const string& name)>        mtllib      = {};
    function<void(const obj_material&)>       material    = {};
    function<void(const obj_camera&)>         camera      = {};
    function<void(const obj_environment&)>    environmnet = {};
};

// Load obj options
struct load_obj_options {
    bool exit_on_error = false;
    bool geometry_only = false;
    bool flip_texcoord = true;
    bool flip_tr       = true;
};

// Load obj scene
bool load_obj(const string& filename, const obj_callbacks& cb,
    const load_obj_options& options = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// ply type
enum struct ply_type { ply_uchar, ply_int, ply_float, ply_int_list };

// ply property
struct ply_property {
    string                name    = "";
    ply_type              type    = ply_type::ply_float;
    vector<float>         scalars = {};
    vector<array<int, 8>> lists   = {};
};

// ply element
struct ply_element {
    string               name       = "";
    int                  count      = 0;
    vector<ply_property> properties = {};
};

// simple ply api data
struct ply_data {
    vector<ply_element> elements = {};
};

// Load ply mesh
bool load_ply(const string& filename, ply_data& ply);

}  // namespace yocto

#endif
