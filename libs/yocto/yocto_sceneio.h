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

#ifndef YOCTO_SCENEIO_H_
#define YOCTO_SCENEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <sstream>
#include <string>
#include <vector>

#include "yocto_ndarray.h"
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
bool is_srgb_filename(const string& filename);

// Loads/saves a float image.
array2d<vec4f> load_image(const string& filename, bool srgb = false);
void           load_image(
              const string& filename, array2d<vec4f>& image, bool srgb = false);
void save_image(
    const string& filename, const array2d<vec4f>& image, bool srgb = false);

// Loads/saves a byte image.
array2d<vec4b> load_imageb(const string& filename, bool srgb = true);
void           load_image(
              const string& filename, array2d<vec4b>& image, bool srgb = true);
void save_image(
    const string& filename, const array2d<vec4b>& image, bool srgb = true);

// Make presets. Supported mostly in IO.
bool           is_srgb_preset(const string& type);
array2d<vec4f> make_image_preset(const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
texture_data load_texture(const string& filename);
void         load_texture(const string& filename, texture_data& texture);
void         save_texture(const string& filename, const texture_data& texture);

// Make presets. Supported mostly in IO.
texture_data make_texture_preset(const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
shape_data load_shape(const string& filename, bool flip_texcoords = true);
void       load_shape(
          const string& filename, shape_data& shape, bool flip_texcoords = true);
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
fvshape_data load_fvshape(const string& filename, bool flip_texcoords = true);
void         load_fvshape(
            const string& filename, fvshape_data& shape, bool flip_texcoords = true);
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
shape_data   make_shape_preset(const string& type);
fvshape_data make_fvshape_preset(const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIV IO
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Using directive
using byte = unsigned char;

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE WATCHER
// -----------------------------------------------------------------------------
namespace yocto {

// File watcher
struct watch_context {
  std::atomic<int>  version   = 0;
  std::future<void> worker    = {};
  vector<string>    filenames = {};
  vector<int64_t>   filetimes = {};
  int64_t           delay     = 500;
  std::atomic<bool> stop      = false;
};

// Initialize file watcher
watch_context make_watch_context(
    const vector<string>& filenames, int delay = 500);
// Start file watcher
void watch_start(watch_context& context);
// Stop file watcher
void watch_stop(watch_context& context);
// Get file version
int get_version(const watch_context& context);

}  // namespace yocto

#endif
