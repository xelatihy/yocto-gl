//
// # Yocto/SceneIO: Scene serialization
//
// Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
// and a custom Json format.
// Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`
// and depends on `stb_image.h`, `stb_image_write.h`, `tinyexr.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
bool load_text(const string& filename, string& str, string& error);
bool save_text(const string& filename, const string& str, string& error);

// Using directive
using byte = unsigned char;

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error);
bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR or LDR based on filename.
bool is_hdr_filename(const string& filename);
bool is_ldr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
bool load_image(const string& filename, color_image& img, string& error);
bool save_image(const string& filename, const color_image& img, string& error);

// Make presets. Supported mostly in IO.
bool make_image_preset(color_image& image, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
bool load_texture(
    const string& filename, scene_texture& texture, string& error);
bool save_texture(
    const string& filename, const scene_texture& texture, string& error);

// Make presets. Supported mostly in IO.
bool make_texture_preset(
    scene_texture& texture, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
bool load_shape(const string& filename, scene_shape& shape, string& error,
    bool flip_texcoords = true);
bool save_shape(const string& filename, const scene_shape& shape, string& error,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
bool load_fvshape(const string& filename, scene_fvshape& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const scene_fvshape& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
bool make_shape_preset(scene_shape& shape, const string& type, string& error);
bool make_fvshape_preset(
    scene_fvshape& shape, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIV IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a subdiv in the supported formats.
bool load_subdiv(const string& filename, scene_subdiv& subdiv, string& error);
bool save_subdiv(
    const string& filename, const scene_subdiv& subdiv, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the supported formats.
// Calls the progress callback, if defined, as we process more data.
bool load_scene(const string& filename, scene_model& scene, string& error,
    bool noparallel = false);
bool save_scene(const string& filename, const scene_model& scene, string& error,
    bool noparallel = false);

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_model& scene, string& error);

// Scene presets used for testing.
bool make_scene_preset(scene_model& scene, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
string normalize_path(const string& filename);

// Get directory name (not including '/').
string path_dirname(const string& filename);

// Get extension (including '.').
string path_extension(const string& filename);

// Get filename without directory.
string path_filename(const string& filename);

// Get filename without directory and extension.
string path_basename(const string& filename);

// Joins paths
string path_join(const string& patha, const string& pathb);
string path_join(const string& patha, const string& pathb, const string& pathc);

// Replaces extensions
string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
bool path_exists(const string& filename);

// Check if a file is a directory
bool path_isdir(const string& filename);

// Check if a file is a file
bool path_isfile(const string& filename);

// List the contents of a directory
vector<string> list_directory(const string& filename);

// Create a directory and all missing parent directories if needed
bool make_directory(const string& dirname, string& error);

// Get the current directory
string path_current();

}  // namespace yocto

#endif
