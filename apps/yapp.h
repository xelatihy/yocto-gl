//
// YAPP: helper scene object to write demo code for YOCTO/GL library.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YAPP_H_
#define _YAPP_H_

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_sym.h"
#include "../yocto/yocto_trace.h"

namespace yapp {

//
// Typedefs for vec/mat types
//
using uint = unsigned int;
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using float4x4 = std::array<std::array<float, 4>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;
using int4 = std::array<int, 4>;

//
// Scene Texture
//
struct texture {
    std::string path;                // path
    int width = 0;                   // width
    int height = 0;                  // height
    int ncomp = 0;                   // number of components
    std::vector<float> hdr;          // if loaded, hdr data
    std::vector<unsigned char> ldr;  // if loaded, ldr data
};

//
// Scene Material
//
struct material {
    // whole material data
    std::string name;  // material name

    // color information
    float3 ke = {0, 0, 0};  // emission color
    float3 kd = {0, 0, 0};  // diffuse color
    float3 ks = {0, 0, 0};  // specular color
    float3 kt = {0, 0, 0};  // transmittance color
    float rs = 0.0001;      // roughness

    // indices in the texture array (-1 if not found)
    texture* ke_txt = nullptr;
    texture* kd_txt = nullptr;
    texture* ks_txt = nullptr;
    texture* kt_txt = nullptr;
    texture* rs_txt = nullptr;
};

//
// Scene geometry
//
struct shape {
    // whole shape data
    std::string name;         // shape name
    material* mat = nullptr;  // material pointer
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame

    // shape elements
    std::vector<int> points;      // points
    std::vector<int2> lines;      // lines
    std::vector<int3> triangles;  // triangles

    // vertex data
    std::vector<float3> pos;       // per-vertex position (3 float)
    std::vector<float3> norm;      // per-vertex normals (3 float)
    std::vector<float2> texcoord;  // per-vertex texcoord (2 float)
    std::vector<float3> color;     // [extension] per-vertex color (3 float)
    std::vector<float> radius;     // [extension] per-vertex radius (1 float)

    // additional vertex data
    std::vector<float3> ke;  // per-vertex emission
    std::vector<float3> kd;  // per-vertex diffuse
    std::vector<float3> ks;  // per-vertex specular
    std::vector<float3> kt;  // per-vertex specular
    std::vector<float> rs;   // per-vertex exponent
};

//
// Scene Camera
//
struct camera {
    std::string name;                                                 // name
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame
    bool ortho = false;                  // ortho camera
    float yfov = std::tan(3.1416f / 3);  // vertical field of view
    float aspect = 16.0f / 9.0f;         // aspect ratio
    float aperture = 0;                  // lens aperture
    float focus = 1;                     // focus distance
};

//
// Envinonment map
//
struct environment {
    std::string name;         // name
    material* mat = nullptr;  // material pointer
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame
};

//
// Asset
//
struct scene {
    std::vector<shape*> shapes;              // shape array
    std::vector<material*> materials;        // material array
    std::vector<texture*> textures;          // texture array
    std::vector<camera*> cameras;            // camera array
    std::vector<environment*> environments;  // environment array

    ~scene();
};

//
// Gets material index
//
int get_material_idx(const scene*& scn, const material*& mat);

int get_etype(const shape& shape);
int get_nelems(const shape& shape);
const int* get_elems(const shape& shape);

//
// Load scene
//
scene* load_scene(const std::string& filename, float scale,
                  bool add_camera = true);

//
// Load scene
//
scene* load_scenes(const std::vector<std::string>& filenames, float scale,
                   bool add_camera = true);

//
// Loads an envmap for the scene
//
void load_envmap(const scene* scn, const std::string& filename, float scale);

//
// Save scene
//
void save_scene(const std::string& filename, const scene* sc);

//
// Make trace blocks
//
std::vector<int4> make_trace_blocks(int w, int h, int bs);

//
// Save image
//
void save_image(const std::string& filename, int width, int height,
                const float4* hdr, float exposure, float gamma,
                bool srgb_output);

//
// Make a BVH
//
ybvh::scene* make_bvh(const scene* scene);

//
// Initialize scene for rendering
//
ytrace::scene* make_trace_scene(const scene* scene,
                                const ybvh::scene* scene_bvh, int camera);

//
// Initialize a rigid body scene
//
ysym::scene* make_rigid_scene(const scene* scene, ybvh::scene*& scene_bvh);

//
// Step one time
//
void simulate_step(scene* scene, ysym::scene* rigid_scene, float dt);

struct params {
    // scene/image
    std::vector<std::string> filenames;
    float scene_scale = 1;

    // render
    std::string imfilename;
    int width = 0, height = 0;
    float exposure = 0;
    float gamma = 1;
    bool srgb = true;
    float4 background = {0, 0, 0, 0};

    // trace
    ytrace::render_params render_params;
    int save_progressive = 0;  // number of frames at which to save
    int block_size = 32;
    int nthreads = 0;

    // rigid
    std::string outfilename;
    float dt = 1 / 60.0f;
    int nframes = 1000;

    // ui
    bool legacy_gl = false;
    bool no_ui = false;

    // shade
    bool wireframe = false, edges = false;
};

//
// Load parameters
//
params* init_params(const std::string& help, int argc, char** argv,
                    bool trace_params, bool sym_params, bool shade_params,
                    bool ui_params);

//
// Logging
//
void set_default_loggers();

// Forward declaration
struct shade_state;

//
// Init OpenGL render
//
shade_state* init_shade_state(const yapp::scene* scn, const params* pars);

//
// Cleanup state
//
void free_shade_state(shade_state* state);

//
// OpenGL render
//
void shade_scene(const yapp::scene* scn, const params* pars,
                 const shade_state* st);

}  // namespace

#endif
