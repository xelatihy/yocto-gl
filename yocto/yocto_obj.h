//
// # Yocto/OBJ: Tiny library for OBJ parsing
//
// Yocto/OBJ is a simple Wavefront OBJ parser that works with callbacks.
// We make no attempt to provide a simple interface for OBJ but just the
// low level parsing code. We support a few extensions such as camera and
// environment map loading.
//
// Error reporting is done through exceptions using the `io_error` exception.
//
// ## Parse an OBJ file
//
// 1. define callbacks in `obj_callback` structure using lambda with capture
//    if desired
// 2. run the parse with `load_obj()`
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_utils.h"

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

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.position == b.position && a.texturecoord == b.texturecoord &&
           a.normal == b.normal;
}

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

    // textures
    obj_texture_info ke_txt;    // emission texture
    obj_texture_info ka_txt;    // ambient texture
    obj_texture_info kd_txt;    // diffuse texture
    obj_texture_info ks_txt;    // specular texture
    obj_texture_info kr_txt;    // reflection texture
    obj_texture_info kt_txt;    // transmission texture
    obj_texture_info ns_txt;    // Phong exponent texture
    obj_texture_info op_txt;    // opacity texture
    obj_texture_info ior_txt;   // ior texture
    obj_texture_info bump_txt;  // bump map

    // pbr values
    bool  has_pbr = false;  // whether pbr values are defined
    float pr      = 0;      // roughness
    float pm      = 0;      // metallic
    float ps      = 0;      // sheen
    float pc      = 0;      // coat
    float pcr     = 0;      // coat roughness

    // textures
    obj_texture_info pr_txt;    // roughness texture
    obj_texture_info pm_txt;    // metallic texture
    obj_texture_info ps_txt;    // sheen texture
    obj_texture_info norm_txt;  // normal map
    obj_texture_info disp_txt;  // displacement map
    obj_texture_info occ_txt;   // occlusion map

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

// Obj procedural object [extension].
struct obj_procedural {
    string  name;                      // name
    frame3f frame = identity_frame3f;  // transform
    string  type;                      // type
    string  material;                  // material
    float   size  = 2;                 // size
    int     level = -1;                // level of subdivision (-1 default)
};

// Obj callbacks
struct obj_callbacks {
    virtual void vert(const vec3f&) {}
    virtual void norm(const vec3f&) {}
    virtual void texcoord(const vec2f&) {}
    virtual void face(const vector<obj_vertex>&) {}
    virtual void line(const vector<obj_vertex>&) {}
    virtual void point(const vector<obj_vertex>&) {}
    virtual void object(const string&) {}
    virtual void group(const string&) {}
    virtual void usemtl(const string&) {}
    virtual void smoothing(const string&) {}
    virtual void mtllib(const string&) {}
    virtual void material(const obj_material&) {}
    virtual void camera(const obj_camera&) {}
    virtual void environmnet(const obj_environment&) {}
    virtual void procedural(const obj_procedural&) {}
};

// Load obj params
struct obj_params {
    bool exit_on_error = false;
    bool geometry_only = false;
    bool flip_texcoord = true;
    bool flip_tr       = true;
};

// Load obj scene
void load_obj(
    const string& filename, obj_callbacks& cb, const obj_params& params = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with unordered_map
template <>
struct hash<yocto::obj_vertex> {
    static constexpr std::hash<int> hasher = std::hash<int>();

    size_t operator()(const yocto::obj_vertex& v) const {
        auto h = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= hasher((&v.position)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

#endif
