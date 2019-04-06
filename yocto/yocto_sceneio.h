//
// # Yocto/SceneIO: Tiny library for Yocto/Scene input and output
//
// Yocto/SceneIO provides loading and saving functionality for scenes
// in Yocto/GL. We support a simple to use YAML format, PLY, OBJ and glTF.
// The YAML serialization is a straight copy of the in-memory scene data.
// To speed up testing, we also support a binary format that is a dump of
// the current scene. This format should not be use for archival though.
//
// Error reporting is done through exceptions using the `io_error` exception.
//
// ## Scene Loading and Saving
//
// 1. load a scene with `load_scene()` and save it with `save_scene()`
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
    bool          skip_meshes               = false;
    bool          obj_split_shapes          = true;
    bool          obj_preserve_face_varying = false;
    atomic<bool>* cancel_flag               = nullptr;
    bool          run_serially              = false;
};
// Scene save options
struct save_scene_options {
    bool          skip_textures = false;
    bool          skip_meshes   = false;
    bool          ply_instances = false;
    atomic<bool>* cancel_flag   = nullptr;
    bool          run_serially  = false;
};

// Load/save a scene in the supported formats.
void load_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
void save_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

}  // namespace yocto

#endif
