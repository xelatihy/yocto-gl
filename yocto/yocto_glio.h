//
// # Yocto/GLIO: Inpuot and output functions for Yocto/GL
//
// Yocto/GLIO provides loading and saving functionality for images and scenes
// in Yocto/GL. For images we support PNG, JPG, TGA, HDR, EXR formats. For
// scene we support a simple to use JSON format, PLY, OBJ and glTF.
// The JSON serialization is a straight copy of the in-memory scene data.
//
// ## Usage
//
// 1. manipulate paths withe path utilities
// 2. load and save text files with `load_text()` and `save_text()`
// 3. load and save binary files with `load_binary()` and `save_binary()`
// 4. load and save images with `load_image()` and `save_image()`
// 5. load a scene with `load_json_scene()` and save it with `save_json_scene()`
// 6. load and save OBJs with `load_obj_scene()` and `save_obj_scene()`
// 7. load and save glTFs with `load_gltf_scene()` and `save_gltf_scene()`
// 6. if desired, the function `load_scene()` and `save_scene()` will either
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

#ifndef _YGLIO_H_
#define _YGLIO_H_

#include "yocto_gl.h"

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Normalize path delimiters.
std::string normalize_path(const std::string& filename);
// Get directory name (not including '/').
std::string get_dirname(const std::string& filename);
// Get extension (not including '.').
std::string get_extension(const std::string& filename);
// Get filename without directory.
std::string get_filename(const std::string& filename);
// Replace extension.
std::string replace_extension(
    const std::string& filename, const std::string& ext);

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Load a text file
std::string load_text(const std::string& filename);
// Save a text file
void save_text(const std::string& filename, const std::string& str);
// Load a binary file
std::vector<byte> load_binary(const std::string& filename);
// Save a binary file
void save_binary(const std::string& filename, const std::vector<byte>& data);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

// Loads/saves a 4 channel image.
image4f load_image(const std::string& filename, float ldr_gamma = 2.2f);
void save_image(
    const std::string& filename, const image4f& img, float ldr_gamma = 2.2f);
image4f load_image_from_memory(
    const byte* data, int data_size, float ldr_gamma = 2.2f);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/save a scene in the supported formats.
std::shared_ptr<scene> load_scene(const std::string& filename,
    bool load_textures = true, bool skip_missing = true);
void save_scene(const std::string& filename, const std::shared_ptr<scene> scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene in the builtin JSON format.
std::shared_ptr<scene> load_json_scene(const std::string& filename,
    bool load_textures = true, bool skip_missing = true);
void save_json_scene(const std::string& filename,
    const std::shared_ptr<scene> scn, bool save_textures = true,
    bool skip_missing = true);

// Load/save a scene from/to OBJ.
std::shared_ptr<scene> load_obj_scene(const std::string& filename,
    bool load_textures = true, bool skip_missing = true,
    bool split_shapes = true);
void save_obj_scene(const std::string& filename,
    const std::shared_ptr<scene> scn, bool save_textures = true,
    bool skip_missing = true);

// Load/save a scene from/to glTF.
std::shared_ptr<scene> load_gltf_scene(const std::string& filename,
    bool load_textures = true, bool skip_missing = true);
void save_gltf_scene(const std::string& filename,
    const std::shared_ptr<scene> scn, bool save_textures = true,
    bool skip_missing = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// MESH IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/Save a mesh
void load_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<vec4f>& color,
    std::vector<float>& radius);
void save_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, const std::vector<vec4f>& color,
    const std::vector<float>& radius, bool ascii = false);

// Load/Save a ply mesh
void load_ply_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<vec4f>& color,
    std::vector<float>& radius);
void save_ply_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, const std::vector<vec4f>& color,
    const std::vector<float>& radius, bool ascii = false);

// Load/Save an OBJ mesh
void load_obj_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, bool flip_texcoord = true);
void save_obj_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, bool flip_texcoord = true);

}  // namespace ygl

#endif
