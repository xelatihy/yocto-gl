//
// # Yocto/SceneIO: Loader and converter for Yocto/Scene
//
// Yocto/SceneIO provides loading and saving functionality for Yocto/Scene
// in a simple to use JSON format. The JSON serialization is a straight copy
// of the in-memory scene data. Textures are stored in standard image formats.
// Meshes are stored in a builtin binary format or as PLY or OBJ.
//
// We also support loading OBJ and glTF formats, but the conversion is best
// effort and likely not robust. For more robust support of these formats,
// please consider Yocto/Obj or Yocto/glTF.
//
// ## Usage
//
// 1. load a scene with `load_json_scene()` and save it with `save_json_scene()`
// 2. convert from and to OBJ with `load_obj_scene()` and `save_obj_scene()`
// 3. convert from and to glTF with `load_gltf_scene()` and `save_gltf_scene()`
// 4. if desired, the function `load_scene()` and `save_scene()` will either
//    load using the internal format or convert on the fly using on the
//    supported conversions
// 5. you can use equivalent functions for meshes
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

#ifndef _YGL_SCENEIO_H_
#define _YGL_SCENEIO_H_

// enable fast parsing
#ifndef YGL_FASTPARSE
#define YGL_FASTPARSE 1
#endif

#include "yocto_math.h"
#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/save a scene in the supported formats.
scene* load_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
void save_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene in the builtin JSON format.
scene* load_json_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
void save_json_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene from/to OBJ.
scene* load_obj_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true, bool split_shapes = true);
void save_obj_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene from/to glTF.
scene* load_gltf_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
void save_gltf_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// MESH IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/Save a mesh
void load_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
void save_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii = false);

// Load/Save a ply mesh
void load_ply_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
void save_ply_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii = false);

// Load/Save an OBJ mesh
void load_obj_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    bool flip_texcoord = true);
void save_obj_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    bool flip_texcoord = true);

}  // namespace ygl

#endif
