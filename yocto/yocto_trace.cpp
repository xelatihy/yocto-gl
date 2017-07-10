//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR YOCTO_TRACE
// -----------------------------------------------------------------------------

#include "yocto_trace.h"

#ifndef YTRACE_NO_BVH
#include "yocto_bvh.h"
#endif

#include <map>

//
// BUG: gltf normalization
// BUG: check recursive pathtraced environment map
// BUG: check __sample_brdf at if (rnl >= wd && rnl < wd + ws) {
// TODO: add envmap sampling and maybe use discreet ditributions
// TODO: check fresnel
//

#include <cassert>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

namespace ytrace {

//
// Camera
//
struct camera {
    ym::frame3f frame = ym::identity_frame3f;  // local-to-world transform
    float yfov = ym::pif / 3;                  // field of view
    float aspect = 1;                          // aspect ratio
    float aperture = 0;                        // lens aperture
    float focus = 1;  // focus plane distance (cannot be zero)
};

//
// Texture
//
struct texture {
    int width = 0;               // width
    int height = 0;              // height
    int ncomp = 0;               // number of components
    const float* hdr = nullptr;  // hdr pixel values;
    const byte* ldr = nullptr;   // ldr pixel values
};

//
// Material type
//
enum struct material_type {
    none,
    emission_only,
    generic,
    gltf_metallic_roughness,
    gltf_specular_glossiness
};

//
// Material
//
struct material {
    // type
    material_type mtype = material_type::none;

    // material values
    ym::vec3f ke = ym::zero3f;  // emission, term
    ym::vec3f kd = ym::zero3f;  // diffuse term
    ym::vec3f ks = ym::zero3f;  // specular term
    ym::vec3f kt = ym::zero3f;  // transmittance term
    float rs = 0.1;             // specular roughness
    float op = 1;               // opacity

    // textures
    texture* ke_txt = nullptr;
    texture* kd_txt = nullptr;
    texture* ks_txt = nullptr;
    texture* kt_txt = nullptr;
    texture* rs_txt = nullptr;
    texture* op_txt = nullptr;
    texture* norm_txt = nullptr;
    texture* occ_txt = nullptr;

    // material flags
    bool use_phong = false;  // whether to use phong
};

//
// Shape
//
struct shape {
    // element data [only one enabled at any given time]
    int nelems = 0;                        // number of elements
    const ym::vec1i* points = nullptr;     // elem data
    const ym::vec2i* lines = nullptr;      // elem data
    const ym::vec3i* triangles = nullptr;  // elem data

    // vertex data
    int nverts = 0;                       // number of vertices
    const ym::vec3f* pos = nullptr;       // vertex data
    const ym::vec3f* norm = nullptr;      // vertex data
    const ym::vec2f* texcoord = nullptr;  // vertex data
    const ym::vec4f* color = nullptr;     // vertex data
    const ym::vec1f* radius = nullptr;    // vertex data
    const ym::vec4f* tangsp = nullptr;    // vertex data

    // per-vertex material
    const ym::vec3f* ke = nullptr;  // vertex data
    const ym::vec3f* kd = nullptr;  // vertex data
    const ym::vec3f* ks = nullptr;  // vertex data
    const ym::vec3f* kt = nullptr;  // vertex data
    const ym::vec1f* rs = nullptr;  // vertex data

    // sampling data
    std::vector<float> cdf;  // for shape, cdf of shape elements for sampling
    float area = 0;          // for shape, shape area
};

//
// Instance
//
struct instance {
    ym::frame3f frame = ym::identity_frame3f;  // local-to-world rigid transform
    material* mat = nullptr;                   // material
    shape* shp = nullptr;                      // shape
};

//
// Environment
//
struct environment {
    ym::frame3f frame = ym::identity_frame3f;  // local-to-world rigid transform
    ym::vec3f ke = ym::zero3f;                 // emission
    texture* ke_txt = nullptr;                 // emission texture
};

//
// Light (either shape or environment).
// This is only used internally and should not be created.
//
struct light {
    instance* ist = nullptr;     // instance
    environment* env = nullptr;  // environment
};

//
// Scene
//
struct scene {
    // intersection callbaks
    void* intersect_ctx = nullptr;                 // ray intersection context
    intersect_first_cb intersect_first = nullptr;  // ray intersection callback
    intersect_any_cb intersect_any = nullptr;      // ray hit callback
#ifndef YTRACE_NO_BVH
    ybvh::scene* intersect_bvh = nullptr;  // intersect internal bvh
#endif

    // scene data
    std::vector<camera*> cameras;            // camera
    std::vector<environment*> environments;  // env
    std::vector<instance*> instances;        // instances
    std::vector<shape*> shapes;              // shapes
    std::vector<material*> materials;        // materials
    std::vector<texture*> textures;          // textures

    // default material
    material* default_material = nullptr;

    // logging callback
    void* logging_ctx = nullptr;           // logging callback
    logging_msg_cb logging_msg = nullptr;  // logging message

    // destructor
    ~scene();

    // [private] light sources
    std::vector<light*> lights;        // lights [private]
    bool shadow_transmission = false;  // wheter to test transmission
};

//
// Scene support
//
scene::~scene() {
    for (auto shp : shapes)
        if (shp) delete shp;
    for (auto ist : instances)
        if (ist) delete ist;
    for (auto mat : materials)
        if (mat) delete mat;
    for (auto txt : textures)
        if (txt) delete txt;
    for (auto cam : cameras)
        if (cam) delete cam;
    for (auto env : environments)
        if (env) delete env;
    for (auto light : lights)
        if (light) delete light;
#ifndef YTRACE_NO_BVH
    if (intersect_bvh) ybvh::free_scene(intersect_bvh);
#endif
}

//
// Logging shortcut
//
static inline void _log(const scene* scn, int level, const char* msg, ...) {
    if (!scn->logging_msg) return;
    va_list args;
    va_start(args, msg);
    scn->logging_msg(level, "yocto_trace", msg, args);
    va_end(args);
}

//
// Public API. See above.
//
scene* make_scene() {
    auto scn = new scene();
    scn->default_material = new material();
    scn->default_material->mtype = material_type::generic;
    scn->default_material->kd = {0.2f, 0.2f, 0.2f};
    return scn;
}

///
// Public API. See above.
///
void free_scene(scene* scn) {
    if (scn) delete scn;
}

//
// Public API. See above.
//
void set_camera(scene* scn, int cid, const ym::frame3f& frame, float yfov,
    float aspect, float aperture, float focus) {
    scn->cameras[cid]->frame = frame;
    scn->cameras[cid]->yfov = yfov;
    scn->cameras[cid]->aspect = aspect;
    scn->cameras[cid]->aperture = aperture;
    scn->cameras[cid]->focus = focus;
}

//
// Public API. See above.
//
int add_camera(scene* scn, const ym::frame3f& frame, float yfov, float aspect,
    float aperture, float focus) {
    scn->cameras.push_back(new camera());
    set_camera(scn, (int)scn->cameras.size() - 1, frame, yfov, aspect);
    return (int)scn->cameras.size() - 1;
}

//
// Public API. See above.
//
void set_texture(
    scene* scn, int tid, int width, int height, int ncomp, const float* hdr) {
    scn->textures[tid]->width = width;
    scn->textures[tid]->height = height;
    scn->textures[tid]->ncomp = ncomp;
    scn->textures[tid]->hdr = hdr;
    scn->textures[tid]->ldr = nullptr;
}

//
// Public API. See above.
//
void set_texture(
    scene* scn, int tid, int width, int height, int ncomp, const byte* ldr) {
    scn->textures[tid]->width = width;
    scn->textures[tid]->height = height;
    scn->textures[tid]->ncomp = ncomp;
    scn->textures[tid]->hdr = nullptr;
    scn->textures[tid]->ldr = ldr;
}

//
// Public API. See above.
//
int add_texture(
    scene* scn, int width, int height, int ncomp, const float* hdr) {
    scn->textures.push_back(new texture());
    set_texture(scn, (int)scn->textures.size() - 1, width, height, ncomp, hdr);
    return (int)scn->textures.size() - 1;
}

//
// Public API. See above.
//
int add_texture(scene* scn, int width, int height, int ncomp, const byte* ldr) {
    scn->textures.push_back(new texture());
    set_texture(scn, (int)scn->textures.size() - 1, width, height, ncomp, ldr);
    return (int)scn->textures.size() - 1;
}

//
// Public API. See above.
//
void set_material_generic(scene* scn, int mid, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, const ym::vec3f& kt, float rs,
    float op, int ke_txt, int kd_txt, int ks_txt, int kt_txt, int rs_txt,
    int op_txt, int norm_txt, int occ_txt, bool use_phong) {
    scn->materials[mid]->mtype = material_type::generic;
    scn->materials[mid]->ke = ke;
    scn->materials[mid]->kd = kd;
    scn->materials[mid]->ks = ks;
    scn->materials[mid]->kt = kt;
    scn->materials[mid]->rs = rs;
    scn->materials[mid]->op = op;
    scn->materials[mid]->ke_txt =
        (ke_txt >= 0) ? scn->textures[ke_txt] : nullptr;
    scn->materials[mid]->kd_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
    scn->materials[mid]->kt_txt =
        (kt_txt >= 0) ? scn->textures[kt_txt] : nullptr;
    scn->materials[mid]->rs_txt =
        (rs_txt >= 0) ? scn->textures[rs_txt] : nullptr;
    scn->materials[mid]->op_txt =
        (op_txt >= 0) ? scn->textures[rs_txt] : nullptr;
    scn->materials[mid]->norm_txt =
        (norm_txt >= 0) ? scn->textures[norm_txt] : nullptr;
    scn->materials[mid]->use_phong = use_phong;
}

//
// Public API. See above.
//
int add_material_generic(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int kt_txt, int rs_txt, int op_txt, int norm_txt,
    int occ_txt, bool use_phong) {
    scn->materials.push_back(new material());
    set_material_generic(scn, (int)scn->materials.size() - 1, ke, kd, ks, kt,
        rs, op, ke_txt, kd_txt, ks_txt, kt_txt, rs_txt, op_txt, norm_txt,
        occ_txt, use_phong);
    return (int)scn->materials.size() - 1;
}

//
// Public API. See above.
//
void set_material(scene* scn, int mid, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, int ke_txt, int kd_txt,
    int ks_txt, int kt_txt, int rs_txt, int norm_txt, bool use_phong) {
    set_material_generic(scn, mid, ke, kd, ks, kt, rs, 1, ke_txt, kd_txt,
        ks_txt, kt_txt, rs_txt, -1, norm_txt, -1, use_phong);
}

//
// Public API. See above.
//
int add_material(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, int ke_txt, int kd_txt,
    int ks_txt, int kt_txt, int rs_txt, int norm_txt, bool use_phong) {
    return add_material_generic(scn, ke, kd, ks, kt, rs, 1, ke_txt, kd_txt,
        ks_txt, kt_txt, rs_txt, norm_txt, -1, use_phong);
}

//
// Public API. See above.
//
void set_material_gltf_metallic_roughness(scene* scn, int mid,
    const ym::vec3f& ke, const ym::vec3f& kd, float ks, float rs, float op,
    int ke_txt, int kd_txt, int ks_txt, int norm_txt, int occ_txt,
    bool use_phong) {
    scn->materials[mid]->mtype = material_type::gltf_metallic_roughness;
    scn->materials[mid]->ke = ke;
    scn->materials[mid]->kd = kd;
    scn->materials[mid]->ks = {ks, ks, ks};
    scn->materials[mid]->rs = rs;
    scn->materials[mid]->op = op;
    scn->materials[mid]->ke_txt =
        (ke_txt >= 0) ? scn->textures[ke_txt] : nullptr;
    scn->materials[mid]->kd_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
    scn->materials[mid]->norm_txt =
        (norm_txt >= 0) ? scn->textures[norm_txt] : nullptr;
    scn->materials[mid]->occ_txt =
        (occ_txt >= 0) ? scn->textures[occ_txt] : nullptr;
    scn->materials[mid]->use_phong = use_phong;
}

//
// Public API. See above.
//
int add_material_gltf_metallic_roughness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, float ks, float rs, float op, int ke_txt, int kd_txt,
    int ks_txt, int norm_txt, int occ_txt, bool use_phong) {
    scn->materials.push_back(new material());
    set_material_gltf_metallic_roughness(scn, (int)scn->materials.size() - 1,
        ke, kd, ks, rs, op, ke_txt, kd_txt, ks_txt, norm_txt, occ_txt,
        use_phong);
    return (int)scn->materials.size() - 1;
}

//
// Public API. See above.
//
void set_material_gltf_specular_glossiness(scene* scn, int mid,
    const ym::vec3f& ke, const ym::vec3f& kd, const ym::vec3f& ks, float rs,
    float op, int ke_txt, int kd_txt, int ks_txt, int norm_txt, int occ_txt,
    bool use_phong) {
    scn->materials[mid]->mtype = material_type::gltf_specular_glossiness;
    scn->materials[mid]->ke = ke;
    scn->materials[mid]->kd = kd;
    scn->materials[mid]->ks = ks;
    scn->materials[mid]->rs = rs;
    scn->materials[mid]->op = op;
    scn->materials[mid]->ke_txt =
        (ke_txt >= 0) ? scn->textures[ke_txt] : nullptr;
    scn->materials[mid]->kd_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
    scn->materials[mid]->norm_txt =
        (norm_txt >= 0) ? scn->textures[norm_txt] : nullptr;
    scn->materials[mid]->occ_txt =
        (occ_txt >= 0) ? scn->textures[occ_txt] : nullptr;
    scn->materials[mid]->use_phong = use_phong;
}

//
// Public API. See above.
//
int add_material_gltf_specular_glossiness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int norm_txt, int occ_txt, bool use_phong) {
    scn->materials.push_back(new material());
    set_material_gltf_specular_glossiness(scn, (int)scn->materials.size() - 1,
        ke, kd, ks, rs, op, ke_txt, kd_txt, ks_txt, norm_txt, occ_txt,
        use_phong);
    return (int)scn->materials.size() - 1;
}

//
// Public API. See above.
//
void set_material_emission_only(scene* scn, int mid, const ym::vec3f& ke,
    int ke_txt, int norm_txt, int occ_txt) {
    scn->materials[mid]->mtype = material_type::emission_only;
    scn->materials[mid]->ke = ke;
    scn->materials[mid]->ke_txt =
        (ke_txt >= 0) ? scn->textures[ke_txt] : nullptr;
    scn->materials[mid]->norm_txt =
        (norm_txt >= 0) ? scn->textures[norm_txt] : nullptr;
    scn->materials[mid]->occ_txt =
        (occ_txt >= 0) ? scn->textures[occ_txt] : nullptr;
}

//
// Public API. See above.
//
int add_material_emission_only(
    scene* scn, const ym::vec3f& ke, int ke_txt, int norm_txt, int occ_txt) {
    scn->materials.push_back(new material());
    set_material_emission_only(
        scn, (int)scn->materials.size() - 1, ke, ke_txt, norm_txt, occ_txt);
    return (int)scn->materials.size() - 1;
}

//
// Public API. See above.
//
void set_environment(scene* scn, int eid, const ym::frame3f& frame,
    const ym::vec3f& ke, int txt_id) {
    scn->environments[eid]->frame = frame;
    scn->environments[eid]->ke = ke;
    scn->environments[eid]->ke_txt =
        (txt_id >= 0) ? scn->textures[txt_id] : nullptr;
}

//
// Public API. See above.
//
int add_environment(
    scene* scn, const ym::frame3f& frame, const ym::vec3f& ke, int txt_id) {
    scn->environments.push_back(new environment());
    set_environment(scn, (int)scn->environments.size() - 1, frame, ke, txt_id);
    return (int)scn->environments.size() - 1;
}

//
// Public API. See above.
//
void set_instance(
    scene* scn, int iid, const ym::frame3f& frame, int sid, int mid) {
    scn->instances[iid]->frame = frame;
    scn->instances[iid]->shp = scn->shapes[sid];
    scn->instances[iid]->mat =
        (mid < 0) ? scn->default_material : scn->materials[mid];
}

//
// Public API. See above.
//
int add_instance(scene* scn, const ym::frame3f& frame, int sid, int mid) {
    scn->instances.push_back(new instance());
    set_instance(scn, (int)scn->instances.size() - 1, frame, sid, mid);
    return (int)scn->instances.size() - 1;
}

//
// Public API. See above.
//
void set_triangle_shape(scene* scn, int sid, int ntriangles,
    const ym::vec3i* triangles, int nverts, const ym::vec3f* pos,
    const ym::vec3f* norm, const ym::vec2f* texcoord, const ym::vec4f* color,
    const ym::vec4f* tangsp) {
    scn->shapes[sid]->nelems = ntriangles;
    scn->shapes[sid]->points = nullptr;
    scn->shapes[sid]->lines = nullptr;
    scn->shapes[sid]->triangles = triangles;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = pos;
    scn->shapes[sid]->norm = norm;
    scn->shapes[sid]->texcoord = texcoord;
    scn->shapes[sid]->color = color;
    scn->shapes[sid]->tangsp = tangsp;
    scn->shapes[sid]->radius = nullptr;
}

//
// Public API. See above.
//
int add_triangle_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord, const ym::vec4f* color,
    const ym::vec4f* tangsp) {
    scn->shapes.push_back(new shape());
    set_triangle_shape(scn, (int)scn->shapes.size() - 1, ntriangles, triangles,
        nverts, pos, norm, texcoord, color, tangsp);
    return (int)scn->shapes.size() - 1;
}

//
// Public API. See above.
//
void set_point_shape(scene* scn, int sid, int npoints, const int* points,
    int nverts, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord, const ym::vec4f* color, const float* radius) {
    scn->shapes[sid]->nelems = npoints;
    scn->shapes[sid]->points = (const ym::vec1i*)points;
    scn->shapes[sid]->lines = nullptr;
    scn->shapes[sid]->triangles = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = pos;
    scn->shapes[sid]->norm = norm;
    scn->shapes[sid]->texcoord = texcoord;
    scn->shapes[sid]->color = color;
    scn->shapes[sid]->radius = (const ym::vec1f*)radius;
    scn->shapes[sid]->tangsp = nullptr;
}

//
// Public API. See above.
//
int add_point_shape(scene* scn, int npoints, const int* points, int nverts,
    const ym::vec3f* pos, const ym::vec3f* norm, const ym::vec2f* texcoord,
    const ym::vec4f* color, const float* radius) {
    scn->shapes.push_back(new shape());
    set_point_shape(scn, (int)scn->shapes.size() - 1, npoints, points, nverts,
        pos, norm, texcoord, color, radius);
    return (int)scn->shapes.size() - 1;
}

//
// Public API. See above.
//
void set_line_shape(scene* scn, int sid, int nlines, const ym::vec2i* lines,
    int nverts, const ym::vec3f* pos, const ym::vec3f* tang,
    const ym::vec2f* texcoord, const ym::vec4f* color, const float* radius) {
    scn->shapes[sid]->nelems = nlines;
    scn->shapes[sid]->points = nullptr;
    scn->shapes[sid]->lines = lines;
    scn->shapes[sid]->triangles = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = pos;
    scn->shapes[sid]->norm = tang;
    scn->shapes[sid]->texcoord = texcoord;
    scn->shapes[sid]->color = color;
    scn->shapes[sid]->radius = (const ym::vec1f*)radius;
    scn->shapes[sid]->tangsp = nullptr;
}

//
// Public API. See above.
//
int add_line_shape(scene* scn, int nlines, const ym::vec2i* lines, int nverts,
    const ym::vec3f* pos, const ym::vec3f* tang, const ym::vec2f* texcoord,
    const ym::vec4f* color, const float* radius) {
    scn->shapes.push_back(new shape());
    set_line_shape(scn, (int)scn->shapes.size() - 1, nlines, lines, nverts, pos,
        tang, texcoord, color, radius);
    return (int)scn->shapes.size() - 1;
}

//
// Sets per-vertex material properties.
//
void set_vert_material(scene* scn, int sid, const ym::vec3f* ke,
    const ym::vec3f* kd, const ym::vec3f* ks, const float* rs) {
    scn->shapes[sid]->ke = (const ym::vec3f*)ke;
    scn->shapes[sid]->kd = (const ym::vec3f*)kd;
    scn->shapes[sid]->ks = (const ym::vec3f*)ks;
    scn->shapes[sid]->rs = (const ym::vec1f*)rs;
}

//
// Sets the intersection callbacks
//
void set_intersection_callbacks(scene* scn, void* ctx,
    intersect_first_cb intersect_first, intersect_any_cb intersect_any) {
    scn->intersect_ctx = ctx;
    scn->intersect_first = intersect_first;
    scn->intersect_any = intersect_any;
}

//
// Sets the logging callbacks
//
void set_logging_callbacks(scene* scn, void* ctx, logging_msg_cb logging_msg) {
    scn->logging_ctx = ctx;
    scn->logging_msg = logging_msg;
}

//
// Phong exponent to roughness. Public API, see above.
//
float specular_exponent_to_roughness(float n) { return sqrtf(2 / (n + 2)); }

//
// Specular to fresnel eta. Public API, see above.
//
void specular_fresnel_from_ks(
    const ym::vec3f& ks, ym::vec3f& es, ym::vec3f& esk) {
    es = {(1 + std::sqrt(ks[0])) / (1 - std::sqrt(ks[0])),
        (1 + std::sqrt(ks[1])) / (1 - std::sqrt(ks[1])),
        (1 + std::sqrt(ks[2])) / (1 - std::sqrt(ks[2]))};
    esk = {0, 0, 0};
}

//
// Compute shape element cdf for shape sampling.
//
template <typename T, typename Weight_callback>
static std::vector<float> _compute_weight_cdf(int nelems, const T* elem,
    float* total_weight, const Weight_callback& weight_cb) {
    // prepare return
    auto cdf = std::vector<float>(nelems);

    // compute weights
    *total_weight = 0;
    for (auto i = 0; i < nelems; i++) {
        auto w = weight_cb(elem[i]);
        *total_weight += w;
        cdf[i] = *total_weight;
    }

    // normalize
    for (auto& c : cdf) c /= *total_weight;

    // done
    return cdf;
}

#ifndef YTRACE_NO_BVH
//
// internal intersection adapter
//
static inline intersect_point _internal_intersect_first(
    void* ctx, const ym::ray3f& ray) {
    auto scene_bvh = (ybvh::scene*)ctx;
    auto isec = ybvh::intersect_scene(scene_bvh, ray, false);
    auto ipt = ytrace::intersect_point();
    ipt.dist = isec.dist;
    ipt.iid = isec.iid;
    ipt.sid = isec.sid;
    ipt.eid = isec.eid;
    ipt.euv = {isec.euv[0], isec.euv[1], isec.euv[2]};
    return ipt;
}

//
// internal intersection adapter
//
static inline bool _internal_intersect_any(void* ctx, const ym::ray3f& ray) {
    auto scene_bvh = (ybvh::scene*)ctx;
    return (bool)ybvh::intersect_scene(scene_bvh, ray, true);
}
#endif

//
// Init acceleation using yocto_bvh. Public API, see above.
//
void init_intersection(scene* scn) {
#ifndef YTRACE_NO_BVH
    scn->intersect_bvh = ybvh::make_scene();
    auto shape_map = std::map<shape*, int>();
    for (auto shp : scn->shapes) {
        if (shp->points) {
            shape_map[shp] = ybvh::add_point_shape(scn->intersect_bvh,
                shp->nelems, (int*)shp->points, shp->nverts, shp->pos,
                (float*)shp->radius);
        } else if (shp->lines) {
            shape_map[shp] =
                ybvh::add_line_shape(scn->intersect_bvh, shp->nelems,
                    shp->lines, shp->nverts, shp->pos, (float*)shp->radius);
        } else if (shp->triangles) {
            shape_map[shp] = ybvh::add_triangle_shape(scn->intersect_bvh,
                shp->nelems, shp->triangles, shp->nverts, shp->pos, nullptr);
        } else {
            assert(false);
        }
    }
    for (auto ist : scn->instances) {
        ybvh::add_instance(scn->intersect_bvh, ist->frame, shape_map[ist->shp]);
    }
    ybvh::build_scene_bvh(scn->intersect_bvh);
    set_intersection_callbacks(scn, scn->intersect_bvh,
        _internal_intersect_first, _internal_intersect_any);
#endif
}

//
// Init lights. Public API, see above.
//
void init_lights(scene* scn) {
    // clear old lights
    for (auto lgt : scn->lights) delete lgt;
    scn->lights.clear();
    scn->shadow_transmission = false;
    for (auto shp : scn->shapes) {
        shp->area = 0;
        shp->cdf.clear();
    }

    for (auto ist : scn->instances) {
        if (ist->mat->kt != ym::zero3f) scn->shadow_transmission = true;
        if (ist->mat->ke == ym::zero3f) continue;
        auto lgt = new light();
        lgt->ist = ist;
        auto shp = ist->shp;
        if (shp->cdf.empty()) {
            if (shp->points) {
                shp->cdf = _compute_weight_cdf(shp->nelems, shp->points,
                    &shp->area, [shp](ym::vec1i e) { return 1; });
            } else if (ist->shp->lines) {
                shp->cdf = _compute_weight_cdf(
                    shp->nelems, shp->lines, &shp->area, [shp](ym::vec2i e) {
                        return ym::length(shp->pos[e[1]] - shp->pos[e[0]]);
                    });
            } else if (shp->triangles) {
                shp->cdf = _compute_weight_cdf(shp->nelems, shp->triangles,
                    &shp->area, [shp](ym::vec3i e) {
                        return ym::triangle_area(
                            shp->pos[e[0]], shp->pos[e[1]], shp->pos[e[2]]);
                    });
            }
        }
        scn->lights.push_back(lgt);
    }

    for (auto env : scn->environments) {
        if (env->ke == ym::zero3f) continue;
        auto lgt = new light();
        lgt->env = env;
        scn->lights.push_back(lgt);
    }
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// Random number smp. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
//
struct _sampler {
    ym::rng_pcg32 rng;  // rnumber number state
    int i, j;           // pixel coordinates
    int s, d;           // sample and dimension indices
    int ns;             // number of samples
    rng_type rtype;     // random number type
};

//
// Initialize a smp ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
//
static inline _sampler _make_sampler(
    int i, int j, int s, int ns, rng_type rtype) {
    // we use various hashes to scramble the pixel values
    _sampler smp = {{0, 0}, i, j, s, 0, ns, rtype};
    uint64_t sample_id = ((uint64_t)(i + 1)) << 0 | ((uint64_t)(j + 1)) << 15 |
                         ((uint64_t)(s + 1)) << 30;
    uint64_t initseq = ym::hash_uint64(sample_id);
    uint64_t initstate =
        ym::hash_uint64(sample_id * 3202034522624059733ull + 1ull);
    init(&smp.rng, initstate, initseq);
    return smp;
}

//
// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
//
static inline float _sample_next1f(_sampler* smp) {
    float rn = 0;
    switch (smp->rtype) {
        case rng_type::uniform: {
            rn = next1f(&smp->rng);
        } break;
        case rng_type::stratified: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            rn = (s + next1f(&smp->rng)) / smp->ns;
        } break;
        case rng_type::cmjs: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            rn = (s + ym::hash_randfloat(s, p * 0xa399d265)) / smp->ns;
        } break;
        default: assert(false);
    }

    smp->d += 1;

    // make sure all sampled numbers are below 1
    if (rn >= 1) rn = 1 - FLT_EPSILON;

    return rn;
}

//
// Generates a 2-dimensional sample.
//
// Implementation notes: see above. Note that using deterministic keyed
// permutaton we can use stratified sampling without preallocating samples.
//
static inline ym::vec2f _sample_next2f(_sampler* smp) {
    ym::vec2f rn = {0, 0};
    switch (smp->rtype) {
        case rng_type::uniform: {
            rn[0] = next1f(&smp->rng);
            rn[1] = next1f(&smp->rng);
        } break;
        case rng_type::stratified: {
            uint32_t ns2 = (uint32_t)std::round(std::sqrt(smp->ns));
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            rn[0] = (s % ns2 + next1f(&smp->rng)) / ns2;
            rn[1] = (s / ns2 + next1f(&smp->rng)) / ns2;
        } break;
        case rng_type::cmjs: {
            uint32_t ns2 = (uint32_t)round(sqrt(smp->ns));
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            int sx = ym::hash_permute(s % ns2, ns2, p * 0xa511e9b3);
            int sy = ym::hash_permute(s / ns2, ns2, p * 0x63d83595);
            float jx = ym::hash_randfloat(s, p * 0xa399d265);
            float jy = ym::hash_randfloat(s, p * 0x711ad6a5);
            rn[0] = (s % ns2 + (sy + jx) / ns2) / ns2;
            rn[1] = (s / ns2 + (sx + jy) / ns2) / ns2;
        } break;
        default: assert(false);
    }

    smp->d += 2;

    // make sure all sampled numbers are below 1
    if (rn[0] >= 1) rn[0] = 1 - FLT_EPSILON;
    if (rn[1] >= 1) rn[1] = 1 - FLT_EPSILON;

    return rn;
}

//
// Creates a 1-dimensional sample in [0,num-1]
//
static inline int _sample_next1i(_sampler* smp, int num) {
    return ym::clamp(int(_sample_next1f(smp) * num), 0, num - 1);
}

//
// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
//
struct point {
    // point type -----------------------------------
    enum struct type {
        none = -1,     // invalid
        env = 0,       // environment
        point = 1,     // points
        line = 2,      // lines
        triangle = 3,  // triangle
    };
    type ptype = type::none;  // element type

    // light id -----------------------------
    const instance* ist = nullptr;     // instance id used for MIS
    const environment* env = nullptr;  // env id used for MIS

    // direction ----------------------------
    ym::vec3f wo = ym::zero3f;  // outgoing direction

    // resolved geometry (shape) ------------
    ym::frame3f frame = ym::identity_frame3f;  // local frame

    // shading ------------------------------
    material_type mtype = material_type::none;  // material type
    ym::vec3f ke = ym::zero3f;                  // material values
    ym::vec3f kd = ym::zero3f;                  // material values
    ym::vec3f ks = ym::zero3f;                  // material values
    ym::vec3f kt = ym::zero3f;                  // material values
    float rs = 0;                               // material values
    float op = 0;                               // material values
    bool use_phong = false;                     // material values

    // helpers ------------------------------
    bool emission_only() const {
        if (ptype == type::none || ptype == type::env) return true;
        return kd == ym::zero3f && ks == ym::zero3f;
    }
};

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static ym::ray3f _eval_camera(
    const camera* cam, const ym::vec2f& uv, const ym::vec2f& luv) {
    auto h = 2 * std::tan(cam->yfov / 2);
    auto w = h * cam->aspect;
    auto o = ym::vec3f{luv[0] * cam->aperture, luv[1] * cam->aperture, 0};
    auto q = ym::vec3f{w * cam->focus * (uv[0] - 0.5f),
        h * cam->focus * (uv[1] - 0.5f), -cam->focus};
    return ym::ray3f(transform_point(cam->frame, o),
        transform_direction(cam->frame, ym::normalize(q - o)));
}

//
// Grab a texture value
//
static inline ym::vec4f _lookup_texture(
    const texture* txt, const ym::vec2i& ij, bool srgb = true) {
    if (txt->ldr) {
        auto v = ym::image_lookup(txt->width, txt->height, txt->ncomp, txt->ldr,
            ij[0], ij[1], (unsigned char)255);
        return (srgb) ? ym::srgb_to_linear(v) : ym::byte_to_float(v);
    } else if (txt->hdr) {
        return ym::image_lookup(
            txt->width, txt->height, txt->ncomp, txt->hdr, ij[0], ij[1], 1.0f);
    } else {
        assert(false);
        return {};
    }
}

//
// Wrapper for above function
//
static ym::vec4f _eval_texture_(
    const texture* txt, const ym::vec2f& texcoord, bool srgb = true) {
    assert(txt);
    assert(txt->hdr || txt->ldr);

    // get image width/height
    auto wh = ym::vec2i{txt->width, txt->height};

    // get coordinates normalized for tiling
    auto st = ym::vec2f{std::fmod(texcoord[0], 1.0f) * wh[0],
        std::fmod(texcoord[1], 1.0f) * wh[1]};
    if (st[0] < 0) st[0] += wh[0];
    if (st[1] < 0) st[1] += wh[1];

    // get image coordinates and residuals
    auto ij = ym::clamp(ym::vec2i{(int)st[0], (int)st[1]}, {0, 0}, wh);
    auto uv = st - ym::vec2f{(float)ij[0], (float)ij[1]};

    // get interpolation weights and indices
    ym::vec2i idx[4] = {ij, {ij[0], (ij[1] + 1) % wh[1]},
        {(ij[0] + 1) % wh[0], ij[1]},
        {(ij[0] + 1) % wh[0], (ij[1] + 1) % wh[1]}};
    auto w = ym::vec4f{(1 - uv[0]) * (1 - uv[1]), (1 - uv[0]) * uv[1],
        uv[0] * (1 - uv[1]), uv[0] * uv[1]};

    // handle interpolation
    return (_lookup_texture(txt, idx[0], srgb) * w[0] +
            _lookup_texture(txt, idx[1], srgb) * w[1] +
            _lookup_texture(txt, idx[2], srgb) * w[2] +
            _lookup_texture(txt, idx[3], srgb) * w[3]);
}

//
// Wrapper for above function
//
static ym::vec4f _eval_texture4(
    const texture* txt, const ym::vec2f& texcoord, bool srgb = true) {
    if (!txt) return {1, 1, 1, 1};
    return _eval_texture_(txt, texcoord, srgb);
}

//
// Wrapper for above function
//
static ym::vec3f _eval_texture3(
    const texture* txt, const ym::vec2f& texcoord, bool srgb = true) {
    if (!txt) return {1, 1, 1};
    auto v = _eval_texture_(txt, texcoord, srgb);
    return {v.x, v.y, v.z};
}

//
// Evaluates emission.
//
static ym::vec3f _eval_emission(const point& pt) {
    if (pt.ke == ym::zero3f) return ym::zero3f;
    switch (pt.ptype) {
        case point::type::env: return pt.ke;
        case point::type::point: return pt.ke;
        case point::type::line: return pt.ke;
        case point::type::triangle:
            return (ym::dot(pt.frame[2], pt.wo) > 0) ? pt.ke : ym::zero3f;
        default: {
            assert(false);
            return ym::zero3f;
        }
    }
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static ym::vec3f _eval_fresnel_dielectric(float cosw, const ym::vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1.0f / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = ym::vec3f{1, 1, 1} - sin2 / eta2;
    if (cos2t[0] < 0 || cos2t[1] < 0 || cos2t[2] < 0)
        return ym::vec3f{1, 1, 1};  // tir

    auto t0 = ym::vec3f{
        std::sqrt(cos2t[0]), std::sqrt(cos2t[1]), std::sqrt(cos2t[2])};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs =
        (ym::vec3f{cosw, cosw, cosw} - t1) / (ym::vec3f{cosw, cosw, cosw} + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2.0f;
}

//
// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static ym::vec3f _eval_fresnel_metal(
    float cosw, const ym::vec3f& eta, const ym::vec3f& etak) {
    if (etak == ym::zero3f) return _eval_fresnel_dielectric(cosw, eta);

    cosw = ym::clamp(cosw, (float)-1, (float)1);
    auto cos2 = cosw * cosw;
    auto sin2 = ym::clamp(1 - cos2, (float)0, (float)1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - ym::vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2 = ym::vec3f{std::sqrt(a2plusb2_2[0]),
        std::sqrt(a2plusb2_2[1]), std::sqrt(a2plusb2_2[2])};
    auto t1 = a2plusb2 + ym::vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a = ym::vec3f{std::sqrt(a_2[0]), std::sqrt(a_2[1]), std::sqrt(a_2[2])};
    auto t2 = 2.0f * a * cosw;
    auto rs = (t1 - t2) / (t1 + t2);

    auto t3 = ym::vec3f{cos2, cos2, cos2} * a2plusb2 +
              ym::vec3f{sin2, sin2, sin2} * ym::vec3f{sin2, sin2, sin2};
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2.0f;
}

//
// Schlick approximation of Fresnel term
//
static ym::vec3f _eval_fresnel_schlick(const ym::vec3f& ks, float cosw) {
    return ks + (ym::vec3f{1, 1, 1} - ks) *
                    std::pow(ym::clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

//
// Evaluates the BRDF scaled by the cosine of the incoming direction.
//
// Implementation notes:
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
// http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
// - uses Kajiya-Kay for hair
// - uses a hack for points
//
static ym::vec3f _eval_brdfcos(const point& pt, const ym::vec3f& wi) {
    // exit if not needed
    if (pt.emission_only()) return ym::zero3f;

    // save wo
    auto wo = pt.wo;

    switch (pt.ptype) {
        case point::type::point: {
            // diffuse term (hack for now)
            auto ido = ym::dot(wo, wi);
            return pt.kd * (2 * ido + 1) / (2 * ym::pif);
        } break;
        case point::type::line: {
            // compute wh
            auto wh = ym::normalize(wo + wi);

            // compute dot products
            auto ndo = ym::dot(pt.frame[2], wo), ndi = ym::dot(pt.frame[2], wi),
                 ndh = ym::clamp(ym::dot(wh, pt.frame[2]), (float)0, (float)1);

            // take sines
            auto so = std::sqrt(ym::clamp(1 - ndo * ndo, (float)0, (float)1)),
                 si = std::sqrt(ym::clamp(1 - ndi * ndi, (float)0, (float)1)),
                 sh = std::sqrt(ym::clamp(1 - ndh * ndh, (float)0, (float)1));

            // exit if needed
            if (si <= 0 || so <= 0) return ym::zero3f;

            // diffuse term (Kajiya-Kay)
            auto diff = pt.kd * si / ym::pif;

            // specular term (Kajiya-Kay)
            auto spec = ym::zero3f;
            if (sh > 0 && pt.ks != ym::zero3f) {
                auto ns = 2 / (pt.rs * pt.rs) - 2;
                auto d = (ns + 2) * std::pow(sh, ns) / (2 + ym::pif);
                spec = pt.ks * si * d / (4.0f * si * so);
            }

            // done
            return diff + spec;
        } break;
        case point::type::triangle: {
            // compute wh
            auto wh = ym::normalize(wo + wi);

            // compute dot products
            auto ndo = ym::dot(pt.frame[2], wo), ndi = ym::dot(pt.frame[2], wi),
                 ndh = ym::clamp(ym::dot(wh, pt.frame[2]), (float)0, (float)1);

            // exit if needed
            if (ndi <= 0 || ndo <= 0) return ym::zero3f;

            // diffuse term
            auto diff = pt.kd * ndi / ym::pif;

            // specular term (GGX)
            auto spec = ym::zero3f;
            if (ndh > 0 && pt.ks != ym::zero3f) {
                // microfacet term
                auto dg = 0.0f;

                if (!pt.use_phong) {
                    // evaluate GGX
                    auto alpha2 = pt.rs * pt.rs;
                    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
                    auto d = alpha2 / (ym::pif * di * di);
#ifndef YTRACE_GGX_SMITH
                    auto lambda_o =
                        (-1 + std::sqrt(
                                  1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) /
                        2;
                    auto lambda_i =
                        (-1 + std::sqrt(
                                  1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) /
                        2;
                    auto g = 1 / (1 + lambda_o + lambda_i);
#else
                    auto go =
                        (2 * ndo) /
                        (ndo + std::sqrt(alpha2 + (1 - alpha2) * ndo * ndo));
                    auto gi =
                        (2 * ndi) /
                        (ndi + std::sqrt(alpha2 + (1 - alpha2) * ndi * ndi));
                    auto g = go * gi;
#endif
                    dg = d * g;
                } else {
                    // evaluate Blinn-Phong
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    dg = (ns + 2) * std::pow(ndh, ns) / (2 + ym::pif);
                }

                // handle fresnel
                auto odh = ym::clamp(dot(wo, wh), 0.0f, 1.0f);
                auto ks = _eval_fresnel_schlick(pt.ks, odh);

                // sum up
                spec = ks * ndi * dg / (4 * ndi * ndo);
            }

            // done
            return diff + spec;
        } break;
        default: {
            assert(false);
            return ym::zero3f;
        }
    }
}

//
// Compute the weight for sampling the BRDF
//
static float _weight_brdfcos(const point& pt, const ym::vec3f& wi) {
    // skip if no component
    if (pt.emission_only()) return 0;

    switch (pt.ptype) {
        case point::type::point:
        case point::type::line: return 4 * ym::pif;
        case point::type::triangle: {
            // save wo
            auto wo = pt.wo;

            // compute wh
            auto wh = ym::normalize(wi + wo);

            // compute dot products
            auto ndo = ym::dot(pt.frame[2], wo), ndi = ym::dot(pt.frame[2], wi),
                 ndh = ym::dot(pt.frame[2], wh);

            // check to make sure we are above the surface
            // updated this for refraction
            if (ndo <= 0 || ndi <= 0) return 0;

            // pick from a sum
            auto wn = ym::max_element_val(pt.kd) + ym::max_element_val(pt.ks);
            auto wd = ym::max_element_val(pt.kd) / wn;
            auto ws = ym::max_element_val(pt.ks) / wn;

            // accumulate probability
            auto pdf = 0.0f;

            // diffuse term
            if (wd && ndi > 0) {
                // homepherical cosine probability
                pdf += wd * ndi / ym::pif;
            }

            // specular term (GGX or Phong)
            if (ws && ndi > 0 && ndo > 0 && ndh > 0) {
                if (!pt.use_phong) {
                    // probability proportional to d * ndh
                    auto cos2 = ndh * ndh;
                    auto tan2 = (1 - cos2) / cos2;
                    auto alpha2 = pt.rs * pt.rs;
                    auto d = alpha2 / (ym::pif * cos2 * cos2 * (alpha2 + tan2) *
                                          (alpha2 + tan2));
                    auto hdo = ym::dot(wo, wh);
                    pdf += ws * d * ndh / (4 * hdo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // compute wh
                    auto wh = ym::normalize(wi + wo);
                    auto ndh = ym::dot(pt.frame[2], wh);
                    // homerispherical cosine power probability
                    pdf += ws * powf(ndh, ns) * (ns + 1) / (2 * ym::pif);
                }
            }

            // done
            return (pdf) ? 1 / pdf : 0;
        } break;
        default: {
            assert(false);
            return 0;
        }
    }
}

//
// Picks a direction based on the BRDF
//
static ym::vec3f _sample_brdfcos(
    const point& pt, float rnl, const ym::vec2f& rn) {
    // skip if no component
    if (pt.emission_only()) return ym::zero3f;

    switch (pt.ptype) {
        case point::type::point:
        case point::type::line: {
            // sample wi with uniform spherical distribution
            auto rz = rn[1], rr = sqrtf(1 - rz * rz),
                 rphi = 2 * ym::pif * rn[0];
            auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
            return transform_direction(pt.frame, wi_local);
        } break;
        case point::type::triangle: {
            // save wo
            auto wo = pt.wo;

            // compute cosine
            auto ndo = ym::dot(pt.frame[2], wo);

            // check to make sure we are above the surface
            // update this for refraction
            if (ndo <= 0) return ym::zero3f;

            // pick from a sum
            auto wn = ym::max_element_val(pt.kd) + ym::max_element_val(pt.ks);
            auto wd = ym::max_element_val(pt.kd) / wn;
            // auto ws = ym::max_element_val(pt.ks) / wn;

            // sample according to diffuse
            if (rnl < wd) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn[1]), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * ym::pif * rn[0];
                // set to wi
                auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(pt.frame, wi_local);
            }

            // sample according to specular (GGX or Phong)
            else {
                // if (rnl >= wd && rnl < wd + ws) {
                if (!pt.use_phong) {
                    // sample wh with ggx distribution
                    auto tan2 = pt.rs * pt.rs * rn[1] / (1 - rn[1]);
                    auto rz = std::sqrt(1 / (tan2 + 1)),
                         rr = std::sqrt(1 - rz * rz),
                         rphi = 2 * ym::pif * rn[0];
                    // set to wh
                    auto wh_local =
                        ym::vec3f{rr * std::cos(rphi), rr * std::sin(rphi), rz};
                    auto wh = transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym::normalize(wh * 2.0f * ym::dot(wo, wh) - wo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // sample wh with hemispherical cosine power distribution
                    auto rz = std::pow(rn[1], 1 / (ns + 1)),
                         rr = std::sqrt(1 - rz * rz),
                         rphi = 2 * ym::pif * rn[0];
                    // set to wh
                    auto wh_local =
                        ym::vec3f{rr * std::cos(rphi), rr * std::sin(rphi), rz};
                    auto wh = transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym::normalize(wh * 2.0f * ym::dot(wo, wh) - wo);
                }
            }
        } break;
        default: {
            assert(false);
            return ym::zero3f;
        }
    }
}

//
// Create a point for an environment map. Resolves material with textures.
//
static point _eval_envpoint(const environment* env, const ym::vec3f& wo) {
    // set shape data
    auto pt = point();

    // env
    pt.env = env;

    // env params
    pt.ptype = point::type::env;

    // direction
    pt.wo = wo;

    // maerial
    pt.ke = env->ke;

    // textures
    if (env->ke_txt) {
        auto w = ym::transform_direction(ym::inverse(env->frame), -wo);
        auto theta =
            (std::acos(ym::clamp(w[1], (float)-1, (float)1)) / ym::pif);
        auto phi = std::atan2(w[2], w[0]) / (2 * ym::pif);
        auto texcoord = ym::vec2f{phi, theta};
        pt.ke *= _eval_texture3(env->ke_txt, texcoord);
    }

    // done
    return pt;
}

//
// Interpolate a value over an element
//
template <int N, int M>
static inline ym::vec<float, N> _interpolate_value(
    const ym::vec<float, N>* vals, const ym::vec<int, M>* elems, int eid,
    const ym::vec3f& euv) {
    auto ret = ym::zero_vec<float, N>();
    if (!vals) return ret;
    auto& elem = elems[eid];
    for (auto i = 0; i < M; i++) ret += vals[elem[i]] * euv[i];
    return ret;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static point _eval_shapepoint(
    const instance* ist, int eid, const ym::vec3f& euv, const ym::vec3f& wo) {
    // set shape data
    auto pt = point();

    // instance
    pt.ist = ist;

    // direction
    pt.wo = wo;

    // shortcuts
    auto shp = ist->shp;
    auto mat = ist->mat;

    // compute points and weights
    auto pos = ym::zero3f, norm = ym::zero3f, ke = ym::zero3f, kd = ym::zero3f,
         ks = ym::zero3f, kt = ym::zero3f;
    auto color = ym::zero4f;
    auto texcoord = ym::zero2f;
    auto tangsp = ym::zero4f;
    auto rs = ym::zero1f;
    if (shp->points) {
        pt.ptype = point::type::point;
        pos = _interpolate_value(shp->pos, shp->points, eid, euv);
        norm =
            ym::normalize(_interpolate_value(shp->norm, shp->points, eid, euv));
        texcoord = _interpolate_value(shp->texcoord, shp->points, eid, euv);
        color = _interpolate_value(shp->color, shp->points, eid, euv);
        ke = _interpolate_value(shp->ke, shp->points, eid, euv);
        kd = _interpolate_value(shp->kd, shp->points, eid, euv);
        ks = _interpolate_value(shp->ks, shp->points, eid, euv);
        kt = _interpolate_value(shp->kt, shp->points, eid, euv);
        rs = _interpolate_value(shp->rs, shp->points, eid, euv);
    } else if (shp->lines) {
        pt.ptype = point::type::line;
        pos = _interpolate_value(shp->pos, shp->lines, eid, euv);
        norm =
            ym::normalize(_interpolate_value(shp->norm, shp->lines, eid, euv));
        texcoord = _interpolate_value(shp->texcoord, shp->lines, eid, euv);
        color = _interpolate_value(shp->color, shp->lines, eid, euv);
        ke = _interpolate_value(shp->ke, shp->lines, eid, euv);
        kd = _interpolate_value(shp->kd, shp->lines, eid, euv);
        ks = _interpolate_value(shp->ks, shp->lines, eid, euv);
        kt = _interpolate_value(shp->kt, shp->lines, eid, euv);
        rs = _interpolate_value(shp->rs, shp->lines, eid, euv);
    } else if (shp->triangles) {
        pt.ptype = point::type::triangle;
        pos = _interpolate_value(shp->pos, shp->triangles, eid, euv);
        norm = _interpolate_value(shp->norm, shp->triangles, eid, euv);
        texcoord = _interpolate_value(shp->texcoord, shp->triangles, eid, euv);
        color = _interpolate_value(shp->color, shp->triangles, eid, euv);
        tangsp = _interpolate_value(shp->tangsp, shp->triangles, eid, euv);
        ke = _interpolate_value(shp->ke, shp->triangles, eid, euv);
        kd = _interpolate_value(shp->kd, shp->triangles, eid, euv);
        ks = _interpolate_value(shp->ks, shp->triangles, eid, euv);
        kt = _interpolate_value(shp->kt, shp->triangles, eid, euv);
        rs = _interpolate_value(shp->rs, shp->triangles, eid, euv);
    }

    // handle normal map
    if (shp->texcoord && shp->tangsp && shp->triangles && mat->norm_txt) {
        auto txt = _eval_texture3(mat->norm_txt, texcoord, false) * 2.0f -
                   ym::vec3f{1, 1, 1};
        auto ntxt = ym::normalize(ym::vec3f{txt[0], -txt[1], txt[2]});
        auto frame = ym::make_frame3_fromzx(
            {0, 0, 0}, norm, {tangsp[0], tangsp[1], tangsp[2]});
        frame[1] *= tangsp[3];
        norm = ym::transform_direction(frame, ntxt);
    }

    // creating frame
    pt.frame = ym::make_frame3_fromz(pos, norm);

    // transform to world space
    pt.frame[3] = ym::transform_point(ist->frame, pt.frame[3]);
    pt.frame[2] = ym::transform_direction(ist->frame, pt.frame[2]);
    pt.frame[0] = ym::transform_direction(ist->frame, pt.frame[0]);
    pt.frame[1] = ym::transform_direction(ist->frame, pt.frame[1]);

    // sample material data
    pt.mtype = ist->mat->mtype;
    pt.ke = mat->ke;
    pt.kd = mat->kd;
    pt.ks = mat->ks;
    pt.kt = mat->kt;
    pt.rs = mat->rs;
    pt.op = mat->op;
    pt.use_phong = mat->use_phong;

    // handle surface color
    if (shp->color) {
        auto c = ym::vec3f{color.x, color.y, color.z};
        pt.ke *= c;
        pt.kd *= c;
        pt.ks *= c;
        pt.kt *= c;
    }

    // handle per-vertex material
    if (shp->ke) pt.ke *= kd;
    if (shp->kd) pt.kd *= kd;
    if (shp->ks) pt.ks *= ks;
    if (shp->kt) pt.kt *= kt;
    if (shp->rs) pt.rs *= rs[0];

    // handle textures
    if (shp->texcoord) {
        if (mat->ke_txt) pt.ke *= _eval_texture3(mat->ke_txt, texcoord);
        switch (mat->mtype) {
            case material_type::none: break;
            case material_type::emission_only: break;
            case material_type::gltf_metallic_roughness: {
                if (mat->kd_txt) {
                    auto kd_txt = _eval_texture4(mat->kd_txt, texcoord);
                    pt.kd *= {kd_txt.x, kd_txt.y, kd_txt.z};
                    pt.op *= kd_txt.w;
                }
                if (mat->ks_txt) {
                    auto ks_txt = _eval_texture4(mat->kd_txt, texcoord);
                    pt.ks *= {ks_txt.z, ks_txt.z, ks_txt.z};
                    pt.rs *= ks_txt.y;
                }
            } break;
            case material_type::gltf_specular_glossiness: {
                if (mat->kd_txt) {
                    auto kd_txt = _eval_texture4(mat->kd_txt, texcoord);
                    pt.kd *= {kd_txt.x, kd_txt.y, kd_txt.z};
                    pt.op *= kd_txt.w;
                }
                if (mat->ks_txt) {
                    auto ks_txt = _eval_texture4(mat->kd_txt, texcoord);
                    pt.ks *= {ks_txt.x, ks_txt.y, ks_txt.z};
                    pt.rs *= ks_txt.w;
                }
            } break;
            case material_type::generic: {
                if (mat->kd_txt) pt.kd *= _eval_texture3(mat->kd_txt, texcoord);
                if (mat->ks_txt) pt.ks *= _eval_texture3(mat->ks_txt, texcoord);
                if (mat->kt_txt) pt.kt *= _eval_texture3(mat->kt_txt, texcoord);
                if (mat->rs_txt)
                    pt.rs *= _eval_texture3(mat->rs_txt, texcoord).x;
            } break;
        }
        if (mat->occ_txt) {
            auto occ = _eval_texture3(mat->occ_txt, texcoord);
            pt.ke *= occ;
            pt.kd *= occ;
            pt.ks *= occ;
            pt.kt *= occ;
        }
    }

    // convert between material types
    switch (mat->mtype) {
        case material_type::none: break;
        case material_type::emission_only: {
            pt.kd = {0, 0, 0};
            pt.ks = {0, 0, 0};
            pt.kt = {0, 0, 0};
            pt.op = 1;
        } break;
        case material_type::gltf_metallic_roughness: {
            auto kb = pt.kd;
            auto ks = ym::vec3f{0.04f, 0.04f, 0.04f};
            auto m = pt.ks.x;
            pt.kd = kb * (1 - m);
            pt.ks = kb * m + ks * (1 - m);
            pt.rs = pt.rs * pt.rs;
            pt.kt = {0, 0, 0};
        } break;
        case material_type::gltf_specular_glossiness: {
            pt.rs = (1 - pt.rs) * (1 - pt.rs);
            pt.kt = {0, 0, 0};
        } break;
        case material_type::generic: break;
    }

    return pt;
}

//
// Sample weight for a light point.
//
static float _weight_light(const point& lpt, const point& pt) {
    switch (lpt.ptype) {
        case point::type::env: {
            return 4 * ym::pif;
        } break;
        case point::type::point: {
            auto d = ym::dist(lpt.frame[3], pt.frame[3]);
            return lpt.ist->shp->area / (d * d);
        } break;
        case point::type::line: {
            assert(false);
            return 0;
        } break;
        case point::type::triangle: {
            auto d = ym::dist(lpt.frame[3], pt.frame[3]);
            return lpt.ist->shp->area *
                   std::abs(ym::dot(lpt.frame[2], lpt.wo)) / (d * d);
        } break;
        default: {
            assert(false);
            return 0;
        } break;
    }
}

//
// Picks a point on a light.
//
static point _sample_light(
    const light* lgt, const point& pt, float rne, const ym::vec2f& rn) {
    if (lgt->ist) {
        auto shp = lgt->ist->shp;
        auto eid = (int)(lower_bound(shp->cdf.begin(), shp->cdf.end(), rne) -
                         shp->cdf.begin());
        if (eid > shp->cdf.size() - 1) eid = (int)shp->cdf.size() - 1;

        auto euv = ym::zero3f;
        if (shp->triangles) {
            euv = {std::sqrt(rn[0]) * (1 - rn[1]), 1 - std::sqrt(rn[0]),
                rn[1] * std::sqrt(rn[0])};
        } else if (shp->lines) {
            euv = {1 - rn[0], rn[0], 0};
        } else if (shp->points) {
            euv = {1, 0, 0};
        } else
            assert(false);
        auto lpt = _eval_shapepoint(lgt->ist, eid, euv, ym::zero3f);
        lpt.wo = ym::normalize(pt.frame[3] - lpt.frame[3]);
        return lpt;
    } else if (lgt->env) {
        auto z = -1 + 2 * rn[1];
        auto rr = std::sqrt(ym::clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * ym::pif * rn[0];
        auto wo = ym::vec3f{std::cos(phi) * rr, z, std::sin(phi) * rr};
        auto lpt = _eval_envpoint(lgt->env, wo);
        return lpt;
    } else {
        assert(false);
        return {};
    }
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f _offset_ray(
    const point& pt, const ym::vec3f& w, const render_params& params) {
    if (dot(w, pt.frame[2]) > 0) {
        return ym::ray3f(
            pt.frame[3] + pt.frame[2] * params.ray_eps, w, params.ray_eps);
    } else {
        return ym::ray3f(
            pt.frame[3] - pt.frame[2] * params.ray_eps, w, params.ray_eps);
    }
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f _offset_ray(
    const point& pt, const point& pt2, const render_params& params) {
    auto ray_dist = (pt2.ptype != point::type::env) ?
                        ym::dist(pt.frame[3], pt2.frame[3]) :
                        FLT_MAX;
    return ym::ray3f(pt.frame[3] + pt.frame[2] * params.ray_eps, -pt2.wo,
        params.ray_eps, ray_dist - 2 * params.ray_eps);
}

//
// Intersects a ray with the scn and return the point (or env point).
//
static point _intersect_scene(const scene* scn, const ym::ray3f& ray) {
    auto isec = scn->intersect_first(scn->intersect_ctx, ray);
    if (isec) {
        return _eval_shapepoint(
            scn->instances[isec.iid], isec.eid, isec.euv, -ray.d);
    } else if (!scn->environments.empty()) {
        return _eval_envpoint(scn->environments[0], -ray.d);
    } else {
        return {};
    }
}

//
// Intersects a scene and offsets the ray
//
static point _intersect_scene(const scene* scn, const point& pt,
    const ym::vec3f& w, const render_params& params) {
    return _intersect_scene(scn, _offset_ray(pt, w, params));
}

//
// Intersects a scene and offsets the ray
//
static inline point _intersect_scene(const scene* scn, const point& pt,
    const point& lpt, const render_params& params) {
    return _intersect_scene(scn, _offset_ray(pt, lpt, params));
}

//
// Transparecy
//
static ym::vec3f _eval_transparency(const point& pt) { return pt.kt; }

//
// Test occlusion
//
static ym::vec3f _eval_transmission(const scene* scn, const point& pt,
    const point& lpt, const render_params& params) {
    if (scn->shadow_transmission) {
        auto cpt = pt;
        auto weight = ym::vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            cpt = _intersect_scene(scn, cpt, lpt, params);
            if (cpt.ptype == point::type::none || cpt.ptype == point::type::env)
                break;
            weight *= _eval_transparency(cpt);
            if (weight == ym::zero3f) break;
        }
        return weight;
    } else {
        auto shadow_ray = _offset_ray(pt, lpt, params);
        return (scn->intersect_any(scn->intersect_ctx, shadow_ray)) ?
                   ym::zero3f :
                   ym::vec3f{1, 1, 1};
    }
}

//
// Mis weight
//
static inline float _weight_mis(float w0, float w1) {
    if (!w0 || !w1) return 1;
    return (1 / w0) / (1 / w0 + 1 / w1);
}

//
// Recursive path tracing.
//
static ym::vec4f _shade_pathtrace(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none ||
        (pt.ptype == point::type::env && params.envmap_invisible))
        return ym::zero4f;

    // emission
    auto l = _eval_emission(pt);
    if (pt.emission_only() || scn->lights.empty()) return {l[0], l[1], l[2], 1};

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // handle transparency
        auto kt = _eval_transparency(pt);
        if (kt != ym::zero3f) {
            auto tprob = ym::max_element_val(kt);
            if (_sample_next1f(smp) < tprob) {
                weight *= kt;
                pt = _intersect_scene(scn, pt, -pt.wo, params);
                emission = true;
                continue;
            }
        }

        // emission
        if (emission) l += weight * _eval_emission(pt);

        // direct – light
        auto lgt = scn->lights[_sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            _sample_light(lgt, pt, _sample_next1f(smp), _sample_next2f(smp));
        auto lld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo) *
                   _weight_light(lpt, pt) * (float)scn->lights.size();
        if (lld != ym::zero3f) {
            l += weight * lld * _eval_transmission(scn, pt, lpt, params);
        }

        // direct – brdf
        auto bpt = _intersect_scene(scn, pt,
            _sample_brdfcos(pt, _sample_next1f(smp), _sample_next2f(smp)),
            params);
        auto bld = _eval_emission(bpt) * _eval_brdfcos(pt, -bpt.wo) *
                   _weight_brdfcos(pt, -bpt.wo);
        if (bld != ym::zero3f) {
            l += weight * bld *
                 _weight_mis(
                     _weight_brdfcos(pt, -bpt.wo), _weight_light(bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;
        if (bpt.emission_only()) break;

        // continue path
        weight *= _eval_brdfcos(pt, -bpt.wo) * _weight_brdfcos(pt, -bpt.wo);
        if (weight == ym::zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rho = pt.kd + pt.ks;
            auto rrprob =
                1.0f -
                std::min(std::max(std::max(rho[0], rho[1]), rho[2]), 0.95f);
            if (_sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        emission = false;
    }

    return {l[0], l[1], l[2], 1};
}

//
// Recursive path tracing.
//
static ym::vec4f _shade_pathtrace_std(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none ||
        (pt.ptype == point::type::env && params.envmap_invisible))
        return ym::zero4f;

    // emission
    auto l = _eval_emission(pt);
    if (pt.emission_only()) return {l[0], l[1], l[2], 1};

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // handle transparency
        auto kt = _eval_transparency(pt);
        if (kt != ym::zero3f) {
            auto tprob = ym::max_element_val(kt);
            if (_sample_next1f(smp) < tprob) {
                weight *= kt;
                pt = _intersect_scene(scn, pt, -pt.wo, params);
                emission = true;
                continue;
            }
        }

        // emission
        if (emission) l += weight * _eval_emission(pt);
        // direct
        auto lgt = scn->lights[_sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            _sample_light(lgt, pt, _sample_next1f(smp), _sample_next2f(smp));
        auto ld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo) *
                  _weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != ym::zero3f) {
            auto shadow_ray = _offset_ray(pt, lpt, params);
            if (!scn->intersect_any(scn->intersect_ctx, shadow_ray))
                l += weight * ld;
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rho = pt.kd + pt.ks;
            auto rrprob =
                1.0f -
                std::min(std::max(std::max(rho[0], rho[1]), rho[2]), 0.95f);
            if (_sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        {
            auto wi =
                _sample_brdfcos(pt, _sample_next1f(smp), _sample_next2f(smp));
            weight *= _eval_brdfcos(pt, wi) * _weight_brdfcos(pt, wi);
            if (weight == ym::zero3f) break;

            pt = _intersect_scene(scn, pt, wi, params);
            emission = false;
            if (pt.emission_only()) break;
        }
    }

    return {l[0], l[1], l[2], 1};
}

//
// Recursive path tracing.
//
static ym::vec4f _shade_pathtrace_hack(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none ||
        (pt.ptype == point::type::env && params.envmap_invisible))
        return ym::zero4f;

    // emission
    auto l = _eval_emission(pt);
    if (pt.emission_only()) return {l[0], l[1], l[2], 1};

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // transmission
        if (pt.kt != ym::zero3f) {
            auto wi = -pt.wo;
            auto ior = 1.4f;
            if (dot(pt.wo, pt.frame[2]) > 0) {
                auto n = pt.frame[2], w = pt.wo;
                auto eta = 1 / ior;
                auto cosi = dot(n, w);
                auto k = 1 - eta * eta * (1 - cosi * cosi);
                assert(k > 0);
                wi = -eta * w + n * (eta * cosi - std::sqrt(k));
                wi = normalize(wi);
            } else {
                auto n = -pt.frame[2], w = pt.wo;
                auto eta = ior;
                auto cosi = dot(n, w);
                auto k = 1 - eta * eta * (1 - cosi * cosi);
                if (k <= 0) break;
                wi = -eta * w + n * (eta * cosi - std::sqrt(k));
                wi = normalize(wi);
            }
            weight *= pt.kt;
            pt = _intersect_scene(scn, pt, wi, params);
            continue;
        }

        // direct
        auto lgt = scn->lights[_sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            _sample_light(lgt, pt, _sample_next1f(smp), _sample_next2f(smp));
        auto ld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo) *
                  _weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != ym::zero3f) {
            auto shadow_ray = _offset_ray(pt, lpt, params);
            if (!scn->intersect_any(scn->intersect_ctx, shadow_ray))
                l += weight * ld;
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rho = pt.kd + pt.ks;
            auto rrprob =
                1.0f -
                std::min(std::max(std::max(rho[0], rho[1]), rho[2]), 0.95f);
            if (_sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        {
            auto wi =
                _sample_brdfcos(pt, _sample_next1f(smp), _sample_next2f(smp));
            weight *= _eval_brdfcos(pt, wi) * _weight_brdfcos(pt, wi);
            if (weight == ym::zero3f) break;

            pt = _intersect_scene(scn, pt, wi, params);
            if (pt.emission_only()) break;
        }
    }

    return {l[0], l[1], l[2], 1};
}

//
// Direct illumination.
//
static ym::vec4f _shade_direct(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none ||
        (pt.ptype == point::type::env && params.envmap_invisible))
        return ym::zero4f;

    // emission
    auto l = _eval_emission(pt);
    if (pt.emission_only()) return {l[0], l[1], l[2], 1};

    // ambient
    l += (ym::vec3f)params.amb * pt.kd;

    // direct
    for (auto& lgt : scn->lights) {
        auto lpt =
            _sample_light(lgt, pt, _sample_next1f(smp), _sample_next2f(smp));
        auto ld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo) *
                  _weight_light(lpt, pt);
        if (ld == ym::zero3f) continue;
        auto shadow_ray = _offset_ray(pt, lpt, params);
        if (scn->intersect_any(scn->intersect_ctx, shadow_ray)) continue;
        l += ld;
    }

    // done
    return {l[0], l[1], l[2], 1};
}

//
// Eyelight for quick previewing.
//
static ym::vec4f _shade_eyelight(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params) {
    // intersection
    point pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none ||
        (pt.ptype == point::type::env && params.envmap_invisible))
        return ym::zero4f;

    // emission
    auto l = _eval_emission(pt);
    if (pt.emission_only()) return {l[0], l[1], l[2], 1};

    // brdf*light
    l += _eval_brdfcos(pt, pt.wo) * ym::pif;

    return {l[0], l[1], l[2], 1};
}

//
// Shader function callback.
//
using shade_fn = ym::vec4f (*)(const scene* scn, const ym::ray3f& ray,
    _sampler* smp, const render_params& params);

//
// Renders a block of pixels. Public API, see above.
//
void trace_block(const scene* scn, int width, int height, ym::vec4f* img,
    int block_x, int block_y, int block_width, int block_height,
    int samples_min, int samples_max, const render_params& params) {
    auto cam = scn->cameras[params.camera_id];
    shade_fn shade;
    switch (params.stype) {
        case shader_type::eyelight: shade = _shade_eyelight; break;
        case shader_type::direct: shade = _shade_direct; break;
        case shader_type::pathtrace: shade = _shade_pathtrace; break;
        default: assert(false); return;
    }
    for (auto j = block_y; j < block_y + block_height; j++) {
        for (auto i = block_x; i < block_x + block_width; i++) {
            auto lp = ym::zero4f;
            for (auto s = samples_min; s < samples_max; s++) {
                auto smp =
                    _make_sampler(i, j, s, params.nsamples, params.rtype);
                auto rn = _sample_next2f(&smp);
                auto uv =
                    ym::vec2f{(i + rn[0]) / width, 1 - (j + rn[1]) / height};
                auto ray = _eval_camera(cam, uv, _sample_next2f(&smp));
                auto l = shade(scn, ray, &smp, params);
                if (!std::isfinite(l[0]) || !std::isfinite(l[1]) ||
                    !std::isfinite(l[2])) {
                    _log(scn, 2, "NaN detected");
                    continue;
                }
                if (params.pixel_clamp > 0)
                    *(ym::vec3f*)&l =
                        ym::clamplen(*(ym::vec3f*)&l, params.pixel_clamp);
                lp += l;
            }
            if (params.progressive && samples_min > 0) {
                img[j * width + i] =
                    (img[j * width + i] * (float)samples_min + lp) /
                    (float)samples_max;
            } else {
                img[j * width + i] = lp / (float)(samples_max - samples_min);
            }
        }
    }
}

//
// Renders the whole image. Public API, see above.
//
void trace_image(const scene* scn, int width, int height, ym::vec4f* pixels,
    const render_params& params) {
    trace_block(scn, width, height, pixels, 0, 0, width, height, 0,
        params.nsamples, params);
}

}  // namespace ytrace

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif
