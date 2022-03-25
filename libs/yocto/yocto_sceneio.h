//
// # Yocto/SceneIO: Scene serialization
//
// Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
// and a custom Json format.
// Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`
// and depends on `stb_image.h`, `stb_image_write.h`, `tinyexr.h`, `cgltf.h`,
// `cgltf_write.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include <sstream>
#include <string>
#include <vector>

#include "yocto_scene.h"
#include "yocto_trace.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IO ERROR
// -----------------------------------------------------------------------------
namespace yocto {

// Result object modeled on std::expected
struct io_error : std::runtime_error {
  using std::runtime_error::runtime_error;
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
bool load_image(const string& filename, image_data& img, string& error);
bool save_image(const string& filename, const image_data& img, string& error);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
image_data load_image(const string& filename);
void       load_image(const string& filename, image_data& image);
void       save_image(const string& filename, const image_data& image);

// Make presets. Supported mostly in IO.
image_data make_image_preset(const string& type);

// Make presets. Supported mostly in IO.
bool make_image_preset(
    const string& filename, image_data& image, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
bool load_texture(const string& filename, texture_data& texture, string& error);
bool save_texture(
    const string& filename, const texture_data& texture, string& error);

// Load/save a texture in the supported formats.
texture_data load_texture(const string& filename);
void         load_texture(const string& filename, texture_data& texture);
void         save_texture(const string& filename, const texture_data& texture);

// Make presets. Supported mostly in IO.
texture_data make_texture_preset(const string& type);

// Make presets. Supported mostly in IO.
bool make_texture_preset(
    const string& filname, texture_data& texture, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a shape
shape_data load_shape(const string& filename, bool flip_texcoords = true);
void       load_shape(
          const string& filename, shape_data& shape, bool flip_texcoords = true);
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
fvshape_data load_fvshape(const string& filename, bool flip_texcoords = true);
void         load_fvshape(
            const string& filename, fvshape_data& shape, bool flip_texcoords = true);
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
shape_data   make_shape_preset(const string& type);
fvshape_data make_fvshape_preset(const string& type);

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
bool load_subdiv(const string& filename, subdiv_data& subdiv, string& error);
bool save_subdiv(
    const string& filename, const subdiv_data& subdiv, string& error);

// Load/save a subdiv in the supported formats.
subdiv_data load_subdiv(const string& filename);
void        load_subdiv(const string& filename, subdiv_data& subdiv);
void        save_subdiv(const string& filename, const subdiv_data& subdiv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the supported formats.
bool load_scene(const string& filename, scene_data& scene, string& error,
    bool noparallel = false);
bool save_scene(const string& filename, const scene_data& scene, string& error,
    bool noparallel = false);

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_data& scene, string& error);

// Add environment
bool add_environment(scene_data& scene, const string& filename, string& error);

// Load/save a scene in the supported formats.
scene_data load_scene(const string& filename, bool noparallel = false);
void       load_scene(
          const string& filename, scene_data& scene, bool noparallel = false);
void save_scene(
    const string& filename, const scene_data& scene, bool noparallel = false);

// Add environment
void add_environment(scene_data& scene, const string& filename);

// Make missing scene directories
void make_scene_directories(const string& filename, const scene_data& scene);

// Scene presets used for testing.
scene_data make_scene_preset(const string& type);

// Scene presets used for testing.
bool make_scene_preset(
    const string& filename, scene_data& scene, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Using directive
using byte = unsigned char;

// Load/save a text file
bool load_text(const string& filename, string& str, string& error);
bool save_text(const string& filename, const string& str, string& error);

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error);
bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

// Load/save a text file
string load_text(const string& filename);
void   load_text(const string& filename, string& str);
void   save_text(const string& filename, const string& str);

// Load/save a binary file
vector<byte> load_binary(const string& filename);
void         load_binary(const string& filename, vector<byte>& data);
void         save_binary(const string& filename, const vector<byte>& data);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PARAMETER IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/Save/Update trace params
trace_params load_trace_params(const string& filename);
void         load_trace_params(const string& filename, trace_params& params);
void         update_trace_params(const string& filename, trace_params& params);
void save_trace_params(const string& filename, const trace_params& params);

// Load/Save/Update color grade params
colorgrade_params load_colorgrade_params(const string& filename);
void load_colorgrade_params(const string& filename, colorgrade_params& params);
void update_colorgrade_params(
    const string& filename, colorgrade_params& params);
void save_colorgrade_params(
    const string& filename, const colorgrade_params& params);

// Load/Save/Update trace params
bool load_trace_params(
    const string& filename, trace_params& params, string& error);
bool update_trace_params(
    const string& filename, trace_params& params, string& error);
bool save_trace_params(
    const string& filename, const trace_params& params, string& error);

// Load/Save/Update color grade params
bool load_colorgrade_params(
    const string& filename, colorgrade_params& params, string& error);
bool update_colorgrade_params(
    const string& filename, colorgrade_params& params, string& error);
bool save_colorgrade_params(
    const string& filename, const colorgrade_params& params, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Get directory name (not including /)
string path_dirname(const string& path);

// Get filename without directory and extension.
string path_basename(const string& path);

// Get extension
string path_extension(const string& path);

// Get path normalized
string path_normalized(const string& path);

// Check if a path exists.
bool path_exists(const string& path);

// Replace the extension of a file
string replace_extension(const string& path, const string& extension);

// Create a directory and all missing parent directories if needed
void make_directory(const string& path);
bool make_directory(const string& path, string& error);

}  // namespace yocto

#endif
