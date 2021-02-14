//
// # Yocto/SceneIO: Scene serialization
//
// Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
// and a custom Json format.
// Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`,
// and depends on `cgltf.h`.
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

#include <functional>
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
using std::function;
using std::pair;
using std::string;
using std::vector;

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

// Load/save a subdiv
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
bool make_shape_preset(shape_data& shape, const string& type, string& error);
bool make_fvshape_preset(
    fvshape_data& shape, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Progress callback called when loading.
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Load/save a scene in the supported formats.
// Calls the progress callback, if defined, as we process more data.
bool load_scene(const string& filename, scene_scene& scene, string& error,
    const progress_callback& progress_cb = {}, bool noparallel = false);
bool save_scene(const string& filename, const scene_scene& scene, string& error,
    const progress_callback& progress_cb = {}, bool noparallel = false);

// Scene presets used for testing.
bool make_scene_preset(scene_scene& scene, const string& type, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ELEMENT IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
bool load_texture(
    const string& filename, scene_texture& texture, string& error);
bool save_texture(
    const string& filename, const scene_texture& texture, string& error);

// Load/save a shape in the supported formats.
bool load_shape(const string& filename, scene_shape& shape, string& error);
bool save_shape(
    const string& filename, const scene_shape& shape, string& error);

// Load/save a subdiv in the supported formats.
bool load_subdiv(const string& filename, scene_subdiv& subdiv, string& error);
bool save_subdiv(
    const string& filename, const scene_subdiv& subdiv, string& error);

}  // namespace yocto

#endif
