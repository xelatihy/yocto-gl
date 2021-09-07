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

#include <string>

#include "yocto_scene.h"

#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IO RESULT
// -----------------------------------------------------------------------------
namespace yocto {

// Result object modeled on std::expected
struct io_status {
  string   error = "";
  explicit operator bool() const { return error.empty(); }
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR or LDR based on filename.
bool is_hdr_filename(const string& filename);
bool is_ldr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
pair<io_status, image_data> load_image(const string& filename);
io_status load_image(const string& filename, image_data& image);
io_status save_image(const string& filename, const image_data& image);

// Make presets. Supported mostly in IO.
image_data make_image_preset(const string& type);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
bool load_image(const string& filename, image_data& img, string& error);
bool save_image(const string& filename, const image_data& img, string& error);

// Make presets. Supported mostly in IO.
bool make_image_preset(
    const string& filename, image_data& image, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
pair<io_status, texture_data> load_texture(const string& filename);
io_status load_texture(const string& filename, texture_data& texture);
io_status save_texture(const string& filename, const texture_data& texture);

// Make presets. Supported mostly in IO.
texture_data make_texture_preset(const string& type);

// Load/save a texture in the supported formats.
bool load_texture(const string& filename, texture_data& texture, string& error);
bool save_texture(
    const string& filename, const texture_data& texture, string& error);

// Make presets. Supported mostly in IO.
bool make_texture_preset(
    const string& filname, texture_data& texture, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
pair<io_status, shape_data> load_shape(
    const string& filename, bool flip_texcoords = true);
io_status load_shape(
    const string& filename, shape_data& shape, bool flip_texcoords = true);
io_status save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
pair<io_status, fvshape_data> load_fvshape(
    const string& filename, bool flip_texcoords = true);
io_status load_fvshape(
    const string& filename, fvshape_data& shape, bool flip_texcoords = true);
io_status save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
shape_data   make_shape_preset(const string& type);
fvshape_data make_fvshape_preset(const string& type);

// Load/save a shape
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
bool make_shape_preset(const string& filname, shape_data& shape, string& error);
bool make_fvshape_preset(
    const string& filname, fvshape_data& shape, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIV IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a subdiv in the supported formats.
pair<io_status, subdiv_data> load_subdiv(const string& filename);
io_status load_subdiv(const string& filename, subdiv_data& subdiv);
io_status save_subdiv(const string& filename, const subdiv_data& subdiv);

// Load/save a subdiv in the supported formats.
bool load_subdiv(const string& filename, subdiv_data& subdiv, string& error);
bool save_subdiv(
    const string& filename, const subdiv_data& subdiv, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the supported formats.
pair<io_status, scene_data> load_scene(
    const string& filename, bool noparallel = false);
io_status load_scene(
    const string& filename, scene_data& scene, bool noparallel = false);
io_status save_scene(
    const string& filename, const scene_data& scene, bool noparallel = false);

// Make missing scene directories
io_status make_scene_directories(
    const string& filename, const scene_data& scene);

// Scene presets used for testing.
scene_data make_scene_preset(const string& type);

// Add environment
io_status add_environment(scene_data& scene, const string& filename);

// Load/save a scene in the supported formats.
bool load_scene(const string& filename, scene_data& scene, string& error,
    bool noparallel = false);
bool save_scene(const string& filename, const scene_data& scene, string& error,
    bool noparallel = false);

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_data& scene, string& error);

// Scene presets used for testing.
bool make_scene_preset(
    const string& filename, scene_data& scene, string& error);

// Add environment
bool add_environment(scene_data& scene, const string& filename, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Using directive
using byte = unsigned char;

// Load/save a text file
pair<io_status, string> load_text(const string& filename);
io_status               load_text(const string& filename, string& str);
io_status               save_text(const string& filename, const string& str);

// Load/save a binary file
pair<io_status, vector<byte>> load_binary(const string& filename);
io_status load_binary(const string& filename, vector<byte>& data);
io_status save_binary(const string& filename, const vector<byte>& data);

// Load/save a text file
bool load_text(const string& filename, string& str, string& error);
bool save_text(const string& filename, const string& str, string& error);

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error);
bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

// Opens a file with utf8 filename
FILE* fopen_utf8(const string& filename, const string& mode);

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

// Get the current directory
string path_current();

// Create a directory and all missing parent directories if needed
bool make_directory(const string& dirname, string& error);

// List the contents of a directory
bool list_directory(
    const string& dirname, vector<string>& entries, string& error);

// Create a directory and all missing parent directories if needed
io_status make_directory(const string& dirname);

// List the contents of a directory
pair<io_status, vector<string>> list_directory(const string& dirname);

}  // namespace yocto

#endif
