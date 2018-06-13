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

#include "../yocto/yocto_gl.h"
#include "../yocto/yocto_glutils.h"

#include <unordered_map>
#include <typeinfo>

namespace ygl {

// Scene selection.
struct scene_selection {
    scene_selection()
        : ptr(nullptr), tinfo(nullptr) {}
    template <typename T>
    scene_selection(const std::shared_ptr<T>& val) : ptr(val), tinfo(&typeid(T)) {}

    template <typename T>
    std::shared_ptr<T> as() const {
        return (&typeid(T) == tinfo) ? std::static_pointer_cast<T>(ptr) : nullptr;
    }

    std::shared_ptr<void> ptr = nullptr;                    // selected pointer
    const std::type_info* tinfo = nullptr;  // type info
};

// Update OpenGL scene data.
void update_gldata(const std::shared_ptr<texture>& txt);
void update_gldata(const std::shared_ptr<shape>& shp);
void update_gldata(const std::shared_ptr<scene>& scn);
void clear_gldata(const std::shared_ptr<texture>& txt);
void clear_gldata(const std::shared_ptr<shape>& shp);
void clear_gldata(const std::shared_ptr<scene>& scn);

// Draw scene with stdsurface program.
void draw_glscene(const std::shared_ptr<scene>& scn, const std::shared_ptr<camera>& cam,
    const glsurface_program& prog, 
    const vec2i& viewport_size, const std::shared_ptr<void>& highlighted, bool eyelight,
    bool wireframe = false, bool edges = false, 
    float exposure = 0, float gamma = 2.2f);

// Handle camera navigation and scene selection
bool handle_glcamera_turntable(const std::shared_ptr<glwindow>& win,  const std::shared_ptr<camera>& cam, float selected_dist);
float handle_glcamera_turntable_dist(const std::shared_ptr<glwindow>& win,  const std::shared_ptr<camera>& cam, const std::shared_ptr<scene>& scn, const scene_selection& sel);
bool handle_glcamera_fps(const std::shared_ptr<glwindow>& win,  const std::shared_ptr<camera>& cam);
bool handle_glscene_selection(const std::shared_ptr<glwindow>& win,  const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, const vec2i& imsize, const frame2f& imframe,
    scene_selection& sel);

// Draws widgets for a camera. Used for quickly making demos.
bool draw_glwidgets_camera_inspector(
    const std::shared_ptr<glwindow>& win, const std::string& lbl, const std::shared_ptr<camera>& cam);

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_glwidgets_scene_tree(const std::shared_ptr<glwindow>& win,  const std::string& lbl, const std::shared_ptr<scene>& scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list, int height = 240,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_glwidgets_scene_inspector(const std::shared_ptr<glwindow>& win,  const std::string& lbl, 
    const std::shared_ptr<scene>& scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list, int height = 240,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

}  // namespace ygl
