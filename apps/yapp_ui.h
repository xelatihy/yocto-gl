//
// Extension to Yocto for application GUIs. Put it ygl namespace for brevity.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include "../yocto/yocto_glutils.h"
#include "../yocto/yocto_scene.h"

#include <unordered_map>

namespace ygl {

// Scene selection.
struct scene_selection {
    scene_selection(std::nullptr_t* val = nullptr)
        : ptr(nullptr), tinfo(nullptr) {}
    template <typename T>
    scene_selection(T* val) : ptr(val), tinfo(&typeid(T)) {}

    template <typename T>
    T* as() {
        return (&typeid(T) == tinfo) ? (T*)ptr : nullptr;
    }

    void* ptr = nullptr;                    // selected pointer
    const std::type_info* tinfo = nullptr;  // type info
};

// Update OpenGL scene data.
void update_gldata(texture* txt);
void update_gldata(shape* shp);
void update_gldata(scene* scn);
void clear_gldata(texture* txt);
void clear_gldata(shape* shp);
void clear_gldata(scene* scn);

// Draw scene with stdsurface program.
void draw_glscene(const scene* scn, const camera* cam,
    const glsurface_program& prog, 
    const vec2i& viewport_size, const void* highlighted, bool eyelight,
    bool wireframe = false, bool edges = false, bool cutout = false,
    float exposure = 0, float gamma = 2.2f, bool cull_backface = false);

// Handle camera navigation and scene selection
bool handle_glcamera_navigation(
    glwindow* win, camera* cam, bool navigation_fps);
bool handle_glscene_selection(glwindow* win, const scene* scn,
    const camera* cam, int res, const vec2f& offset, float zoom,
    scene_selection& sel);

// Draws widgets for a camera. Used for quickly making demos.
bool draw_imgui_camera_inspector(
    glwindow* win, const std::string& lbl, camera* cam);

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_imgui_scene_tree(glwindow* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_imgui_scene_inspector(glwindow* win, const std::string& lbl,
    scene* scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

}  // namespace ygl
