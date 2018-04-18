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
    scene_selection() : ptr(nullptr), tinfo(nullptr) {}
    template <typename T>
    scene_selection(T* val) : ptr(val), tinfo(&typeid(T)) {}

    template <typename T>
    T* as() {
        return (&typeid(T) == tinfo) ? (T*)ptr : nullptr;
    }

    void* ptr = nullptr;                    // selected pointer
    const std::type_info* tinfo = nullptr;  // type info
};

// Vertex buffers for scene drawing. Members are not part of the public API.
struct glshape {
    glvertex_buffer pos;        // position
    glvertex_buffer norm;       // normals
    glvertex_buffer texcoord;   // texcoord
    glvertex_buffer texcoord1;  // texcoord
    glvertex_buffer tangsp;     // tangent space
    glvertex_buffer color;      // color
    glvertex_buffer points;     // point elements
    glvertex_buffer lines;      // line elements
    glvertex_buffer triangles;  // triangle elements
    glvertex_buffer quads;      // quad elements as tris
    glvertex_buffer beziers;    // bezier elements as l.
    glvertex_buffer edges;      // edge elements
};

// Initialize gl lights.
gllights make_gllights(const scene* scn);

// Update scene textures on the GPU.
void update_gltexture(const texture* txt, gltexture& gtxt);

// Update scene textures on the GPU.
inline std::unordered_map<texture*, gltexture> make_gltextures(
    const scene* scn) {
    auto gtextures = std::unordered_map<texture*, gltexture>();
    for (auto txt : scn->textures) update_gltexture(txt, gtextures[txt]);
    return gtextures;
}

// Clear OpenGL state.
inline void clear_gltextures(std::unordered_map<texture*, gltexture>& txts) {
    for (auto& kv : txts) clear_gltexture(kv.second);
    txts.clear();
}

// Update scene shapes on the GPU.
void update_glshape(const shape* shp, glshape& gshp);

// Clear OpenGL state.
void clear_glshape(glshape& gshp);

// Update scene shapes on the GPU.
inline std::unordered_map<shape*, glshape> make_glshapes(const scene* scn) {
    auto gshapes = std::unordered_map<shape*, glshape>();
    for (auto shp : scn->shapes) update_glshape(shp, gshapes[shp]);
    return gshapes;
}

// Clear OpenGL state.
inline void clear_glshapes(std::unordered_map<shape*, glshape>& shps) {
    for (auto& kv : shps) clear_glshape(kv.second);
    shps.clear();
}

// Params for stdsurface drawing.
struct glsurface_params {
    int resolution = 512;               // image resolution
    bool wireframe = false;             // wireframe drawing
    bool edges = false;                 // draw edges
    float edge_offset = 0.01f;          // offset for edges
    bool cutout = false;                // draw with binary transparency
    bool eyelight = false;              // camera light mode
    float exposure = 0;                 // exposure
    float gamma = 2.2f;                 // gamma
    vec4f background = {0, 0, 0, 0};    // background color
    vec3f ambient = {0, 0, 0};          // ambient lighting
    vec3f highlight_color = {1, 1, 0};  // highlight color
    vec3f edge_color = {0, 0, 0};       // edge color
    bool double_sided = false;          // double sided rendering
    bool cull_backface = false;         // culling back face
};

// Draw scene with stdsurface program.
void draw_glsurface_scene(const scene* scn, const camera* cam,
    glsurface_program& prog, std::unordered_map<shape*, glshape>& shapes,
    std::unordered_map<texture*, gltexture>& textures, const gllights& lights,
    const vec2i& viewport_size, const void* highlighted,
    const glsurface_params& params);

// Handle camera navigation and scene selection
bool handle_glcamera_navigation(
    glwindow* win, camera* cam, bool navigation_fps);
bool handle_glscene_selection(glwindow* win, const scene* scn,
    const camera* cam, const bvh_tree* bvh, int res,
    const vec2f& offset, float zoom, scene_selection& sel);

// Draws widgets for params.
bool draw_imgui_stdsurface_inspector(
    glwindow* win, const std::string& lbl, glsurface_params& params);

// Draws a widget that can selected the camera.
inline bool draw_imgui_camera_selector(glwindow* win, const std::string& lbl,
    camera*& cam, const scene* scn, camera* view) {
    return draw_imgui_combobox(win, lbl, cam, scn->cameras, true, view);
}

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
