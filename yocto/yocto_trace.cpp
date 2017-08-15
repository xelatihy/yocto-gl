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

#include "yocto_utils.h"

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
#pragma GCC diagnostic ignored "-Wformat-security"
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
    int width = 0;                   // width
    int height = 0;                  // height
    const ym::vec4f* hdr = nullptr;  // hdr pixel values;
    const ym::vec4b* ldr = nullptr;  // ldr pixel values
};

//
// Reflectance type
//
enum struct reflectance_type {
    none,
    matte,
    microfacet,
    gltf_metallic_roughness,
    gltf_specular_glossiness,
    thin_glass,
};

// material reflectance - matte
struct reflectance_matte {
    ym::vec3f kd = ym::zero3f;  // diffuse term
    float op = 1;               // opacity

    texture* kd_txt = nullptr;  // diffuse texture
    texture* op_txt = nullptr;  // opacity texture
};

// material reflectance - microfacet
struct reflectance_microfacet {
    ym::vec3f kd = ym::zero3f;  // diffuse term
    ym::vec3f ks = ym::zero3f;  // specular term
    ym::vec3f kt = ym::zero3f;  // transmission term
    float rs = 0.1;             // specular roughness
    float op = 1;               // opacity
    bool use_phong = false;     // whether to use phong

    texture* kd_txt = nullptr;  // diffuse texture
    texture* ks_txt = nullptr;  // specular texture
    texture* kt_txt = nullptr;  // transmission texture
    texture* rs_txt = nullptr;  // roughness texture
    texture* op_txt = nullptr;  // opacity texture
};

// material reflectance - gltf_metallic_roughness
struct reflectance_gltf_metallic_roughness {
    ym::vec3f kb = ym::zero3f;  // base term
    float km = 0;               // metallic term
    float rs = 0.1;             // specular roughness
    float op = 1;               // opacity

    texture* kb_txt = nullptr;  // base texture
    texture* km_txt = nullptr;  // metallic texture
};

// material reflectance - gltf_specular_glossiness
struct reflectance_gltf_specular_glossiness {
    ym::vec3f kd = ym::zero3f;  // diffuse term
    ym::vec3f ks = ym::zero3f;  // specular term
    float rs = 0.1;             // specular roughness
    float op = 1;               // opacity

    texture* kd_txt = nullptr;  // diffuse texture
    texture* ks_txt = nullptr;  // specular texture
};

// material reflectance - thin glass
struct reflectance_thin_glass {
    ym::vec3f ks = ym::zero3f;  // specular term
    ym::vec3f kt = ym::zero3f;  // transmission term

    texture* ks_txt = nullptr;  // specular texture
    texture* kt_txt = nullptr;  // transmission texture
};

//
// Material
//
struct material {
    // type
    reflectance_type rtype = reflectance_type::none;

    // material emission
    ym::vec3f ke = ym::zero3f;  // emission term
    texture* ke_txt = nullptr;  // emission texture

    // material reflectance - matte
    reflectance_matte matte;
    // material reflectance - microfacet
    reflectance_microfacet microfacet;
    // material reflectance - gltf_metallic_roughness
    reflectance_gltf_metallic_roughness metalrough;
    // material reflectance - gltf_specular_glossiness
    reflectance_gltf_specular_glossiness specgloss;
    // material reflectance - thin glass
    reflectance_thin_glass thin_glass;

    // other textures
    texture* norm_txt = nullptr;  // nromal texture
    texture* occ_txt = nullptr;   // occlusion texture

    // flags
    bool double_sided = false;  // double sided

    // conservative test for opacity
    bool is_opaque() const {
        switch (rtype) {
            case reflectance_type::none: return false;
            case reflectance_type::matte: return matte.op != 1 || matte.op_txt;
            case reflectance_type::microfacet:
                return microfacet.op != 1 || microfacet.op_txt ||
                       microfacet.kt != ym::zero3f || microfacet.kt_txt;
            case reflectance_type::gltf_metallic_roughness:
                return metalrough.op != 1 || metalrough.kb_txt;
            case reflectance_type::gltf_specular_glossiness:
                return specgloss.op != 1 || specgloss.kd_txt;
            case reflectance_type::thin_glass:
                return microfacet.kt != ym::zero3f || microfacet.kt_txt;
        }
    }
};

//
// Shape
//
struct shape {
    // element data [only one enabled at any given time]
    int nelems = 0;                        // number of elements
    const int* points = nullptr;           // elem data
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
    logging_cb log_info = nullptr;   // logging message
    logging_cb log_error = nullptr;  // logging message

    // destructor
    ~scene() {
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

    // [private] light sources
    std::vector<light*> lights;        // lights [private]
    bool shadow_transmission = false;  // wheter to test transmission
};

//
// Public API. See above.
//
scene* make_scene() {
    auto scn = new scene();
    scn->default_material = new material();
    scn->default_material->rtype = reflectance_type::matte;
    scn->default_material->matte.kd = {0.2f, 0.2f, 0.2f};
    return scn;
}

///
// Public API. See above.
///
void free_scene(scene*& scn) {
    if (scn) delete scn;
    scn = nullptr;
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
    scene* scn, int tid, int width, int height, const ym::vec4f* hdr) {
    scn->textures[tid]->width = width;
    scn->textures[tid]->height = height;
    scn->textures[tid]->hdr = hdr;
    scn->textures[tid]->ldr = nullptr;
}

//
// Public API. See above.
//
void set_texture(
    scene* scn, int tid, int width, int height, const ym::vec4b* ldr) {
    scn->textures[tid]->width = width;
    scn->textures[tid]->height = height;
    scn->textures[tid]->hdr = nullptr;
    scn->textures[tid]->ldr = ldr;
}

//
// Public API. See above.
//
int add_texture(scene* scn, int width, int height, const ym::vec4f* hdr) {
    scn->textures.push_back(new texture());
    set_texture(scn, (int)scn->textures.size() - 1, width, height, hdr);
    return (int)scn->textures.size() - 1;
}

//
// Public API. See above.
//
int add_texture(scene* scn, int width, int height, const ym::vec4b* ldr) {
    scn->textures.push_back(new texture());
    set_texture(scn, (int)scn->textures.size() - 1, width, height, ldr);
    return (int)scn->textures.size() - 1;
}

//
// Public API. See above.
//
int add_material(scene* scn) {
    scn->materials.push_back(new material());
    return (int)scn->materials.size() - 1;
}

//
// Public API. See above.
//
void set_material_emission(
    scene* scn, int mid, const ym::vec3f& ke, int ke_txt) {
    scn->materials[mid]->ke = ke;
    scn->materials[mid]->ke_txt =
        (ke_txt >= 0) ? scn->textures[ke_txt] : nullptr;
}

//
// Public API. See above.
//
void set_material_normal(scene* scn, int mid, int norm_txt, float scale) {
    scn->materials[mid]->norm_txt =
        (norm_txt >= 0) ? scn->textures[norm_txt] : nullptr;
}

//
// Public API. See above.
//
void set_material_occlusion(scene* scn, int mid, int occ_txt, float scale) {
    scn->materials[mid]->occ_txt =
        (occ_txt >= 0) ? scn->textures[occ_txt] : nullptr;
}

//
// Public API. See above.
//
void set_material_microfacet(scene* scn, int mid, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, float op, int kd_txt,
    int ks_txt, int kt_txt, int rs_txt, int op_txt, bool use_phong) {
    scn->materials[mid]->rtype = reflectance_type::microfacet;
    scn->materials[mid]->microfacet.kd = kd;
    scn->materials[mid]->microfacet.ks = ks;
    scn->materials[mid]->microfacet.kt = kt;
    scn->materials[mid]->microfacet.rs = rs;
    scn->materials[mid]->microfacet.op = op;
    scn->materials[mid]->microfacet.kd_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->microfacet.ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
    scn->materials[mid]->microfacet.kt_txt =
        (kt_txt >= 0) ? scn->textures[kt_txt] : nullptr;
    scn->materials[mid]->microfacet.rs_txt =
        (rs_txt >= 0) ? scn->textures[rs_txt] : nullptr;
    scn->materials[mid]->microfacet.op_txt =
        (op_txt >= 0) ? scn->textures[rs_txt] : nullptr;
    scn->materials[mid]->microfacet.use_phong = use_phong;
}

//
// Public API. See above.
//
void set_material_gltf_metallic_roughness(scene* scn, int mid,
    const ym::vec3f& kb, float km, float rs, float op, int kd_txt, int km_txt) {
    scn->materials[mid]->rtype = reflectance_type::gltf_metallic_roughness;
    scn->materials[mid]->metalrough.kb = kb;
    scn->materials[mid]->metalrough.km = km;
    scn->materials[mid]->metalrough.rs = rs;
    scn->materials[mid]->metalrough.op = op;
    scn->materials[mid]->metalrough.kb_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->metalrough.km_txt =
        (km_txt >= 0) ? scn->textures[km_txt] : nullptr;
}

//
// Public API. See above.
//
void set_material_gltf_specular_glossiness(scene* scn, int mid,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op, int kd_txt,
    int ks_txt) {
    scn->materials[mid]->rtype = reflectance_type::gltf_specular_glossiness;
    scn->materials[mid]->specgloss.kd = kd;
    scn->materials[mid]->specgloss.ks = ks;
    scn->materials[mid]->specgloss.rs = rs;
    scn->materials[mid]->specgloss.op = op;
    scn->materials[mid]->specgloss.kd_txt =
        (kd_txt >= 0) ? scn->textures[kd_txt] : nullptr;
    scn->materials[mid]->specgloss.ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
}

//
// Public API. See above.
//
void set_material_thin_glass(scene* scn, int mid, const ym::vec3f& ks,
    const ym::vec3f& kt, int ks_txt, int kt_txt) {
    scn->materials[mid]->rtype = reflectance_type::thin_glass;
    scn->materials[mid]->thin_glass.ks = ks;
    scn->materials[mid]->thin_glass.kt = kt;
    scn->materials[mid]->thin_glass.ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
    scn->materials[mid]->thin_glass.ks_txt =
        (ks_txt >= 0) ? scn->textures[ks_txt] : nullptr;
}

//
// Double sided material
//
void set_material_double_sided(scene* scn, int mid, bool double_sided) {
    scn->materials[mid]->double_sided = double_sided;
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
    scn->shapes[sid]->points = points;
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
// Sets the intersection callbacks
//
void set_intersection_callbacks(scene* scn, intersect_first_cb intersect_first,
    intersect_any_cb intersect_any) {
    scn->intersect_first = intersect_first;
    scn->intersect_any = intersect_any;
}

//
// Sets the logging callbacks
//
void set_logging_callbacks(
    scene* scn, logging_cb logging_info, logging_cb logging_error) {
    scn->log_info = logging_info;
    scn->log_error = logging_error;
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
    es = {(1 + std::sqrt(ks.x)) / (1 - std::sqrt(ks.x)),
        (1 + std::sqrt(ks.y)) / (1 - std::sqrt(ks.y)),
        (1 + std::sqrt(ks.z)) / (1 - std::sqrt(ks.z))};
    esk = {0, 0, 0};
}

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
    set_intersection_callbacks(scn,
        [scn](const ym::ray3f& ray) {
            auto isec = ybvh::intersect_scene(scn->intersect_bvh, ray, false);
            auto ipt = ytrace::intersect_point();
            ipt.dist = isec.dist;
            ipt.iid = isec.iid;
            ipt.sid = isec.sid;
            ipt.eid = isec.eid;
            ipt.euv = {isec.euv.x, isec.euv.y, isec.euv.z};
            return ipt;
        },
        [scn](const ym::ray3f& ray) {
            return (bool)ybvh::intersect_scene(scn->intersect_bvh, ray, true);
        });
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
        if (!ist->mat->is_opaque()) scn->shadow_transmission = true;
        if (ist->mat->ke == ym::zero3f) continue;
        auto lgt = new light();
        lgt->ist = ist;
        auto shp = ist->shp;
        if (shp->cdf.empty()) {
            shp->cdf.resize(shp->nelems);
            if (shp->points) {
                shp->cdf.resize(shp->nelems);
                for (auto i = 0; i < shp->nelems; i++) ist->shp->cdf[i] = i + 1;
            } else if (ist->shp->lines) {
                ym::sample_lines_cdf(
                    shp->nelems, shp->lines, shp->pos, shp->cdf.data());
            } else if (shp->triangles) {
                ym::sample_triangles_cdf(
                    shp->nelems, shp->triangles, shp->pos, shp->cdf.data());
            }
            shp->area = ist->shp->cdf.back();
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
struct sampler {
    ym::rng_pcg32 rng;  // rnumber number state
    int i, j;           // pixel coordinates
    int s, d;           // sample and dimension indices
    int ns, ns2;        // number of samples and its square root
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
static inline sampler make_sampler(
    int i, int j, int s, int ns, rng_type rtype) {
    // we use various hashes to scramble the pixel values
    sampler smp = {
        {0, 0}, i, j, s, 0, ns, (int)std::round(std::sqrt((float)ns)), rtype};
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
static inline float sample_next1f(sampler* smp) {
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
static inline ym::vec2f sample_next2f(sampler* smp) {
    ym::vec2f rn = {0, 0};
    switch (smp->rtype) {
        case rng_type::uniform: {
            rn.x = next1f(&smp->rng);
            rn.y = next1f(&smp->rng);
        } break;
        case rng_type::stratified: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            rn.x = (s % smp->ns2 + next1f(&smp->rng)) / smp->ns2;
            rn.y = (s / smp->ns2 + next1f(&smp->rng)) / smp->ns2;
        } break;
        case rng_type::cmjs: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(smp->i + 1)) << 0 |
                                            ((uint64_t)(smp->j + 1)) << 15 |
                                            ((uint64_t)(smp->d + 1)) << 30);
            int s = ym::hash_permute(smp->s, smp->ns, p);
            int sx = ym::hash_permute(s % smp->ns2, smp->ns2, p * 0xa511e9b3);
            int sy = ym::hash_permute(s / smp->ns2, smp->ns2, p * 0x63d83595);
            float jx = ym::hash_randfloat(s, p * 0xa399d265);
            float jy = ym::hash_randfloat(s, p * 0x711ad6a5);
            rn.x = (s % smp->ns2 + (sy + jx) / smp->ns2) / smp->ns2;
            rn.y = (s / smp->ns2 + (sx + jy) / smp->ns2) / smp->ns2;
        } break;
        default: assert(false);
    }

    smp->d += 2;

    // make sure all sampled numbers are below 1
    if (rn.x >= 1) rn.x = 1 - FLT_EPSILON;
    if (rn.y >= 1) rn.y = 1 - FLT_EPSILON;

    return rn;
}

//
// Creates a 1-dimensional sample in [0,num-1]
//
static inline int sample_next1i(sampler* smp, int num) {
    return ym::clamp(int(sample_next1f(smp) * num), 0, num - 1);
}

//
// Brdf type
//
enum struct brdf_type {
    reflection_lambert = 0,
    reflection_ggx = 1,
    transmission_lambert = 10,
    transmission_ggx = 11,
    refraction_ggx = 21,
    kajiya_kay_diff = 30,
    kajiya_kay_spec = 31,
    point_diffuse = 40,
    transparent = 50
};

//
// Brdf lobe
//
struct brdf {
    brdf_type type = brdf_type::reflection_lambert;
    ym::vec3f rho = ym::zero3f;
    float roughness = 1;
};

//
// Emission type
//
enum struct emission_type {
    diffuse = 0,
    point = 1,
    line = 2,
    env = 3,
};

//
// Emission lobe
//
struct emission {
    emission_type type = emission_type::diffuse;
    ym::vec3f ke = ym::zero3f;
};

//
// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
//
struct point {
    // light id -----------------------------
    const instance* ist = nullptr;     // instance id used for MIS
    const environment* env = nullptr;  // env id used for MIS

    // direction ----------------------------
    ym::vec3f wo = ym::zero3f;  // outgoing direction

    // resolved geometry (shape) ------------
    ym::frame3f frame = ym::identity_frame3f;  // local frame

    // shading ------------------------------
    float op = 0;                // material opacity
    int nemissions = 0;          // numner of emission lobes
    emission emissions[8] = {};  // emission lobes
    int nbrdfs = 0;              // numner of brdf lobes
    brdf brdfs[8] = {};          // brdf lobes
    ym::vec3f rho = ym::zero3f;  // material brdf weight
    ym::vec2f texcoord;          // texcoord for debugging

    // helpers ------------------------------
    bool no_reflectance() const { return nbrdfs == 0; }
    bool no_emission() const { return nemissions == 0; }
};

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static ym::ray3f eval_camera(
    const camera* cam, const ym::vec2f& uv, const ym::vec2f& luv) {
    auto h = 2 * std::tan(cam->yfov / 2);
    auto w = h * cam->aspect;
    auto o = ym::vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    auto q = ym::vec3f{w * cam->focus * (uv.x - 0.5f),
        h * cam->focus * (uv.y - 0.5f), -cam->focus};
    return ym::ray3f(transform_point(cam->frame, o),
        transform_direction(cam->frame, ym::normalize(q - o)));
}

//
// Grab a texture value
//
static inline ym::vec4f lookup_texture(
    const texture* txt, const ym::vec2i& ij, bool srgb = true) {
    if (txt->ldr) {
        auto v = txt->ldr[ij.y * txt->width + ij.x];
        return (srgb) ? ym::srgb_to_linear(v) : ym::byte_to_float(v);
    } else if (txt->hdr) {
        return txt->hdr[ij.y * txt->width + ij.x];
    } else {
        assert(false);
        return {};
    }
}

//
// Wrapper for above function
//
static ym::vec4f eval_texture(
    const texture* txt, const ym::vec2f& texcoord, bool srgb = true) {
    if (!txt) return {1, 1, 1, 1};
    assert(txt->hdr || txt->ldr);

    // get image width/height
    auto wh = ym::vec2i{txt->width, txt->height};

    // get coordinates normalized for tiling
    auto st = ym::vec2f{
        std::fmod(texcoord.x, 1.0f) * wh.x, std::fmod(texcoord.y, 1.0f) * wh.y};
    if (st.x < 0) st.x += wh.x;
    if (st.y < 0) st.y += wh.y;

    // get image coordinates and residuals
    auto ij = ym::clamp(ym::vec2i{(int)st.x, (int)st.y}, {0, 0}, wh);
    auto uv = st - ym::vec2f{(float)ij.x, (float)ij.y};

    // get interpolation weights and indices
    ym::vec2i idx[4] = {ij, {ij.x, (ij.y + 1) % wh.y},
        {(ij.x + 1) % wh.x, ij.y}, {(ij.x + 1) % wh.x, (ij.y + 1) % wh.y}};
    auto w = ym::vec4f{(1 - uv.x) * (1 - uv.y), (1 - uv.x) * uv.y,
        uv.x * (1 - uv.y), uv.x * uv.y};

    // handle interpolation
    return (lookup_texture(txt, idx[0], srgb) * w.x +
            lookup_texture(txt, idx[1], srgb) * w.y +
            lookup_texture(txt, idx[2], srgb) * w.z +
            lookup_texture(txt, idx[3], srgb) * w.w);
}

//
// Evaluates emission.
//
static ym::vec3f eval_emission(const point& pt) {
    if (!pt.nemissions) return ym::zero3f;
    auto ke = ym::zero3f;
    for (auto lid = 0; lid < pt.nemissions; lid++) {
        auto emission = pt.emissions[lid];
        switch (emission.type) {
            case emission_type::diffuse:
                ke +=
                    (ym::dot(pt.frame.z, pt.wo) > 0) ? emission.ke : ym::zero3f;
                break;
            case emission_type::point: ke += emission.ke; break;
            case emission_type::line: ke += emission.ke; break;
            case emission_type::env: ke += emission.ke; break;
            default: assert(false); break;
        }
    }
    return ke;
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static ym::vec3f eval_fresnel_dielectric(float cosw, const ym::vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1.0f / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = ym::vec3f{1, 1, 1} - sin2 / eta2;
    if (cos2t.x < 0 || cos2t.y < 0 || cos2t.z < 0)
        return ym::vec3f{1, 1, 1};  // tir

    auto t0 =
        ym::vec3f{std::sqrt(cos2t.x), std::sqrt(cos2t.y), std::sqrt(cos2t.z)};
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
static ym::vec3f eval_fresnel_metal(
    float cosw, const ym::vec3f& eta, const ym::vec3f& etak) {
    if (etak == ym::zero3f) return eval_fresnel_dielectric(cosw, eta);

    cosw = ym::clamp(cosw, (float)-1, (float)1);
    auto cos2 = cosw * cosw;
    auto sin2 = ym::clamp(1 - cos2, (float)0, (float)1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - ym::vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2 = ym::vec3f{std::sqrt(a2plusb2_2.x), std::sqrt(a2plusb2_2.y),
        std::sqrt(a2plusb2_2.z)};
    auto t1 = a2plusb2 + ym::vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a = ym::vec3f{std::sqrt(a_2.x), std::sqrt(a_2.y), std::sqrt(a_2.z)};
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
static ym::vec3f eval_fresnel_schlick(const ym::vec3f& ks, float cosw) {
    return ks + (ym::vec3f{1, 1, 1} - ks) *
                    std::pow(ym::clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

//
// Schlick approximation of Fresnel term weighted by roughness.
// This is a hack, but works better than not doing it.
//
static ym::vec3f eval_fresnel_schlick(
    const ym::vec3f& ks, float cosw, float rs) {
    auto fks = eval_fresnel_schlick(ks, cosw);
    return ym::lerp(ks, fks, rs);
}

//
// Evaluates the GGX distribution and geometric term
//
static inline float eval_ggx(float rs, float ndh, float ndi, float ndo) {
    // evaluate GGX
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (ym::pif * di * di);
#ifndef YTRACE_GGX_SMITH
    auto lambda_o =
        (-1 + std::sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i =
        (-1 + std::sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    auto go = (2 * ndo) / (ndo + std::sqrt(alpha2 + (1 - alpha2) * ndo * ndo));
    auto gi = (2 * ndi) / (ndi + std::sqrt(alpha2 + (1 - alpha2) * ndi * ndi));
    auto g = go * gi;
#endif
    return d * g;
}

//
// Evaluates the GGX pdf
//
static inline float pdf_ggx(float rs, float ndh) {
    auto cos2 = ndh * ndh;
    auto tan2 = (1 - cos2) / cos2;
    auto alpha2 = rs * rs;
    auto d =
        alpha2 / (ym::pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
    return d;
}

//
// Sample the GGX distribution
//
static ym::vec3f sample_ggx(float rs, const ym::vec2f& rn) {
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz = std::sqrt(1 / (tan2 + 1)), rr = std::sqrt(1 - rz * rz),
         rphi = 2 * ym::pif * rn.x;
    // set to wh
    auto wh_local = ym::vec3f{rr * std::cos(rphi), rr * std::sin(rphi), rz};
    return wh_local;
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
static ym::vec3f eval_brdfcos(const point& pt, const ym::vec3f& wi) {
    // exit if not needed
    if (pt.no_reflectance()) return ym::zero3f;

    // save wo
    auto wo = pt.wo;

    // accumulate brdfcos for each lobe
    auto brdfcos = ym::zero3f;
    // keep a weight to account for fresnel term
    auto weight = ym::vec3f{1, 1, 1};
    for (auto lid = pt.nbrdfs - 1; lid >= 0; lid--) {
        auto brdf = pt.brdfs[lid];
        switch (brdf.type) {
            // reflection terms
            case brdf_type::reflection_lambert:
            case brdf_type::reflection_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // exit if needed
                if (ndi <= 0 || ndo <= 0) continue;

                // diffuse term
                if (brdf.type == brdf_type::reflection_lambert) {
                    brdfcos += weight * brdf.rho * ndi / ym::pif;
                }
                // specular term
                else {
                    // compute wh
                    auto wh = ym::normalize(wo + wi);

                    // compute dot products
                    auto ndh =
                        ym::clamp(ym::dot(wh, pt.frame.z), (float)-1, (float)1);

                    // exit if needed
                    if (ndh <= 0) continue;

                    // microfacet term
                    auto dg = eval_ggx(brdf.roughness, ndh, ndi, ndo);

                    // handle fresnel
                    auto odh = ym::clamp(dot(wo, wh), 0.0f, 1.0f);
                    auto ks =
                        eval_fresnel_schlick(brdf.rho, odh, brdf.roughness);

                    // sum up
                    brdfcos += weight * ks * ndi * dg / (4 * ndi * ndo);

                    // update weight
                    weight *= ym::vec3f{1, 1, 1} - ks;
                }
            } break;
            // transmission terms
            case brdf_type::transmission_lambert:
            case brdf_type::transmission_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // exit if needed
                if (ndi >= 0 || ndo <= 0) continue;

                // flip direction
                ndi = -ndi;
                auto wi_ = -wi;

                // diffuse term
                if (brdf.type == brdf_type::reflection_lambert) {
                    brdfcos += weight * brdf.rho * ndi / ym::pif;
                }
                // specular term
                else {
                    // compute wh
                    auto wh = ym::normalize(wo + wi_);

                    // compute dot products
                    auto ndh =
                        ym::clamp(ym::dot(wh, pt.frame.z), (float)-1, (float)1);

                    // exit if needed
                    if (ndh <= 0) continue;

                    // microfacet term
                    auto dg = eval_ggx(brdf.roughness, ndh, ndi, ndo);

                    // handle fresnel
                    auto odh = ym::clamp(dot(wo, wh), 0.0f, 1.0f);
                    auto ks =
                        eval_fresnel_schlick(brdf.rho, odh, brdf.roughness);

                    // sum up
                    brdfcos += weight * ks * ndi * dg / (4 * ndi * ndo);
                }
            } break;
            // refraction term (GGX)
            case brdf_type::refraction_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // HACK
                // eta
                auto eta = 1.4f;

                // HACK
                auto kt = brdf.rho;

                // exit if needed
                if (ndi * ndo >= 0) continue;

                // flip eta if necessary
                if (ndo < 0) eta = 1 / eta;

                // compute wh
                auto wh = ym::normalize(wo + eta * wi);

                // flip from pbrt
                if (dot(pt.frame.z, wh) < 0) wh = -wh;

                // compute dot products
                auto ndh = ym::clamp(
                         ym::dot(wh, pt.frame.z), (float)-1, (float)1),
                     odh = ym::dot(wh, wo), idh = ym::dot(wh, wi);

                // microfacet term
                auto dg = eval_ggx(brdf.roughness, std::fabs(ndh),
                    std::fabs(ndi), std::fabs(ndo));

                // fresnel
                auto f = 1 - eval_fresnel_schlick({0.04, 0.04, 0.04}, ndh).x;

                // sum up
                brdfcos += weight * kt * std::fabs(ndi) *
                           (std::fabs(idh) * std::fabs(odh) * f * dg) /
                           (std::fabs(ndi) * std::fabs(ndo) *
                               (eta * idh + odh) * (eta * idh + odh));

            } break;
            // hair (Kajiya-Kay)
            case brdf_type::kajiya_kay_diff:
            case brdf_type::kajiya_kay_spec: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // take sines
                auto so = std::sqrt(
                         ym::clamp(1 - ndo * ndo, (float)0, (float)1)),
                     si = std::sqrt(
                         ym::clamp(1 - ndi * ndi, (float)0, (float)1));

                // exit if needed
                if (si <= 0 || so <= 0) continue;

                // diffuse term (Kajiya-Kay)
                if (brdf.type == brdf_type::kajiya_kay_diff) {
                    brdfcos += weight * brdf.rho * si / ym::pif;
                }
                // specular term (Kajiya-Kay)
                else {
                    // compute wh
                    auto wh = ym::normalize(wo + wi);

                    // compute dot products
                    auto ndh =
                        ym::clamp(ym::dot(wh, pt.frame.z), (float)0, (float)1);

                    // take sines
                    auto sh =
                        std::sqrt(ym::clamp(1 - ndh * ndh, (float)0, (float)1));

                    // exit if needed
                    if (sh <= 0) continue;

                    // specular
                    auto ns = 2 / (brdf.roughness * brdf.roughness) - 2;
                    auto d = (ns + 2) * std::pow(sh, ns) / (2 + ym::pif);
                    brdfcos += weight * brdf.rho * si * d / (4.0f * si * so);
                }
            } break;
            // points
            case brdf_type::point_diffuse: {
                auto ido = ym::dot(wo, wi);
                brdfcos += weight * brdf.rho * (2 * ido + 1) / (2 * ym::pif);
            } break;
            // transparent
            case brdf_type::transparent: {
                // compute cosines
                auto ido = ym::dot(wo, wi);

                // exit if needed
                if (ido > -0.999f) continue;

                // transparent transmission hack
                brdfcos += weight * brdf.rho;
            } break;
            default: assert(false); break;
        }
    }

    // check
    assert(ym::isfinite(brdfcos));

    // done
    return brdfcos;
}

//
// Compute the weight for sampling the BRDF
//
static float weight_brdfcos(const point& pt, const ym::vec3f& wi) {
    // skip if no component
    if (pt.no_reflectance()) return 0;

    // probability of each lobe
    auto weights = ym::vec<float, 8>();
    auto sum = 0.0f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        weights[lid] = ym::max_element_val(pt.brdfs[lid].rho);
        sum += weights[lid];
    }
    if (!sum) return 0;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) weights[lid] /= sum;

    // save wo
    auto wo = pt.wo;

    // accumulate the probability over all lobes
    auto pdf = 0.0f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        auto brdf = pt.brdfs[lid];

        // sample the lobe
        switch (brdf.type) {
            // reflection term
            case brdf_type::reflection_lambert:
            case brdf_type::reflection_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // check to make sure we are above the surface
                if (ndo <= 0 || ndi <= 0) continue;

                // diffuse term
                if (brdf.type == brdf_type::reflection_lambert) {
                    // hemipherical cosine probability
                    pdf += weights[lid] * ndi / ym::pif;
                }
                // specular term (GGX)
                else {
                    // compute wh
                    auto wh = ym::normalize(wi + wo);

                    // compute dot products
                    auto ndh = ym::dot(pt.frame.z, wh);

                    // check to make sure we are above the surface
                    if (ndh <= 0) continue;

                    // specular term (GGX)
                    // probability proportional to d adjusted by wh projection
                    auto d = pdf_ggx(brdf.roughness, ndh);
                    auto hdo = ym::dot(wo, wh);
                    pdf += weights[lid] * d / (4 * hdo);
                }
            } break;
            // transmission term
            case brdf_type::transmission_lambert:
            case brdf_type::transmission_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // check to make sure we are above the surface
                if (ndo <= 0 || ndi >= 0) continue;

                // flip
                ndi = -ndi;
                auto wi_ = -wi;

                // diffuse term
                if (brdf.type == brdf_type::reflection_lambert) {
                    // hemipherical cosine probability
                    pdf += weights[lid] * ndi / ym::pif;
                }
                // specular term (GGX)
                else {
                    // compute wh
                    auto wh = ym::normalize(wi_ + wo);

                    // compute dot products
                    auto ndh = ym::dot(pt.frame.z, wh);

                    // check to make sure we are above the surface
                    if (ndh <= 0) continue;

                    // specular term (GGX)
                    // probability proportional to d adjusted by wh projection
                    auto d = pdf_ggx(brdf.roughness, ndh);
                    auto hdo = ym::dot(wo, wh);
                    pdf += weights[lid] * d / (4 * hdo);
                }
            } break;
            // refraction term (GGX)
            case brdf_type::refraction_ggx: {
                // compute dot products
                auto ndo = ym::dot(pt.frame.z, wo),
                     ndi = ym::dot(pt.frame.z, wi);

                // exit if needed
                if (ndi * ndo >= 0) continue;

                // HACK
                // eta
                auto eta = 1.4f;

                // flip eta if necessary
                if (ndo < 0) eta = 1 / eta;

                // compute wh
                auto wh = ym::normalize(wo + eta * wi);

                // compute dot products
                auto ndh = ym::clamp(
                         ym::dot(wh, pt.frame.z), (float)-1, (float)1),
                     odh = ym::dot(wh, wo), idh = ym::dot(wh, wi);

                // specular term (GGX)
                // probability proportional to d weighted by change of variable
                auto d = pdf_ggx(brdf.roughness, std::fabs(ndh));
                // pdf += weights[lid] * d * eta * eta * std::fabs(idh) /
                //        ((odh + eta * idh) * (odh + eta * idh));

                static ym::vec2f acc = ym::zero2f;
                static ym::vec2i count = ym::zero2i;

                auto x = d * std::fabs(idh) /
                         ((odh + eta * idh) * (odh + eta * idh));
                auto idx = (ndo < 0) ? 0 : 1;
                acc[idx] += x;
                count[idx] += 1;
                pdf += weights[lid] * d * std::fabs(idh) /
                       ((odh + eta * idh) * (odh + eta * idh));

                // check
                assert(ym::isfinite(pdf));
            } break;
            // diffuse term (Kajiya-Kay)
            case brdf_type::kajiya_kay_diff:
            // specular term (Kajiya-Kay)
            case brdf_type::kajiya_kay_spec:
            // diffuse term point
            case brdf_type::point_diffuse: {
                pdf += weights[lid] * 4 * ym::pif;
            } break;
            // transparent term point
            case brdf_type::transparent: {
                // compute dot products
                auto ido = ym::dot(wi, wo);

                // check to make sure we are resonably close to trasmission
                if (ido > 0.999f) continue;

                // add probability
                pdf += weights[lid];
            } break;
            default: assert(false); break;
        }
    }

    // check for missed pdf
    if (!pdf) return 0;

    // check
    assert(ym::isfinite(pdf));

    // done
    return 1 / pdf;
}

//
// reflected vector
//
static inline ym::vec3f reflect(const ym::vec3f& w, const ym::vec3f& n) {
    return -w + 2 * dot(n, w) * n;
}

//
// refracted vector
//
static inline ym::vec3f refract(
    const ym::vec3f& w, const ym::vec3f& n, float eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * std::max(0.0f, 1 - ym::dot(n, w) * ym::dot(n, w));
    if (k < 0) return ym::zero3f;  // tir
    return -w * eta + (eta * ym::dot(n, w) - std::sqrt(k)) * n;
}

//
// Picks a direction based on the BRDF
//
static ym::vec3f sample_brdfcos(
    const point& pt, float rnl, const ym::vec2f& rn) {
    // skip if no component
    if (pt.no_reflectance()) return ym::zero3f;

    // probability of each lobe
    auto weights = ym::vec<float, 8>();
    auto sum = 0.0f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        weights[lid] = ym::max_element_val(pt.brdfs[lid].rho);
        sum += weights[lid];
    }
    if (!sum) return ym::zero3f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) weights[lid] /= sum;

    // pick a lobe
    auto cdf = ym::vec<float, 8>();
    for (auto lid = 0; lid < pt.nbrdfs; lid++)
        cdf[lid] = weights[lid] + ((lid > 0) ? cdf[lid - 1] : 0.0f);
    auto lid = 0;
    for (lid = 0; lid < 8; lid++)
        if (rnl < cdf[lid]) break;
    lid = ym::clamp(lid, 0, pt.nbrdfs - 1);
    auto brdf = pt.brdfs[lid];

    // value to be returned
    auto wi = ym::zero3f;

    // sample selected lobe
    switch (brdf.type) {
        // reflection term
        case brdf_type::reflection_lambert:
        case brdf_type::reflection_ggx: {
            // save wo
            auto wo = pt.wo;

            // compute cosine
            auto ndo = ym::dot(pt.frame.z, wo);

            // check to make sure we are above the surface
            if (ndo <= 0) return ym::zero3f;

            // sample according to diffuse
            if (brdf.type == brdf_type::reflection_lambert) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * ym::pif * rn.x;
                // set to wi
                auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(pt.frame, wi_local);
            }
            // sample according to specular GGX
            else {
                // sample wh with ggx distribution
                auto wh_local = sample_ggx(brdf.roughness, rn);
                auto wh = transform_direction(pt.frame, wh_local);
                // compute wi
                return ym::normalize(wh * 2.0f * ym::dot(wo, wh) - wo);
            }
        } break;
        // tranbsmission term
        case brdf_type::transmission_lambert:
        case brdf_type::transmission_ggx: {
            // save wo
            auto wo = pt.wo;

            // compute cosine
            auto ndo = ym::dot(pt.frame.z, wo);

            // check to make sure we are above the surface
            if (ndo <= 0) return ym::zero3f;

            // sample according to diffuse
            if (brdf.type == brdf_type::reflection_lambert) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * ym::pif * rn.x;
                // set to wi
                auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return -transform_direction(pt.frame, wi_local);
            }
            // sample according to specular GGX
            else {
                // sample wh with ggx distribution
                auto wh_local = sample_ggx(brdf.roughness, rn);
                auto wh = transform_direction(pt.frame, wh_local);
                // compute wi
                return -ym::normalize(wh * 2.0f * ym::dot(wo, wh) - wo);
            }
        } break;
        // transmission term (GGX)
        case brdf_type::refraction_ggx: {
            // save wo
            auto wo = pt.wo;

            // compute cosine
            auto ndo = ym::dot(pt.frame.z, wo);

            // HACK
            // eta
            auto eta = 1.4f;

            // flip eta if necessary
            if (ndo < 0) eta = 1 / eta;

            // sample according to specular (GGX or Phong)
            // sample wh with ggx distribution
            auto wh_local = sample_ggx(brdf.roughness, rn);
            auto wh = transform_direction(pt.frame, wh_local);

            // wi
            auto e = 1 / eta;
            auto wi = ym::zero3f;
            auto odh = dot(pt.frame.z, wo);
            auto k = 1 - e * e * std::max(0.0f, 1 - odh * odh);
            if (k < 0)
                wi = ym::zero3f;  // tir
            else if (ndo < 0) {
                wi = ym::normalize(-wo * e + (e * odh + std::sqrt(k)) * wh);
                assert(dot(pt.frame.z, wi) * ndo <= 0);

            } else {
                wi = ym::normalize(-wo * e + (e * odh - std::sqrt(k)) * wh);
                assert(dot(pt.frame.z, wi) * ndo <= 0);
            }

            // check
            assert(ym::isfinite(wi));

            // done
            return wi;
        } break;
        // diffuse term (Kajiya-Kay)
        case brdf_type::kajiya_kay_diff:
        // specular term (Kajiya-Kay)
        case brdf_type::kajiya_kay_spec:
        // diffuse term point
        case brdf_type::point_diffuse: {
            // sample wi with uniform spherical distribution
            auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz),
                 rphi = 2 * ym::pif * rn.x;
            auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
            return transform_direction(pt.frame, wi_local);
        } break;
        // transparent term
        case brdf_type::transparent: {
            // continue ray direction
            return -pt.wo;
        } break;
        default: assert(false); break;
    }

    // done
    return wi;
}  // namespace ytrace

//
// Create a point for an environment map. Resolves material with textures.
//
static point eval_envpoint(const environment* env, const ym::vec3f& wo) {
    // set shape data
    auto pt = point();

    // env
    pt.env = env;

    // direction
    pt.wo = wo;

    // maerial
    auto ke = env->ke;
    if (env->ke_txt) {
        auto w = ym::transform_direction(ym::inverse(env->frame), -wo);
        auto theta = (std::acos(ym::clamp(w.y, (float)-1, (float)1)) / ym::pif);
        auto phi = std::atan2(w.z, w.x) / (2 * ym::pif);
        auto texcoord = ym::vec2f{phi, theta};
        ke *= eval_texture(env->ke_txt, texcoord).xyz();
    }

    // create emission lobe
    if (ke != ym::zero3f) {
        pt.emissions[0].type = emission_type::env;
        pt.emissions[0].ke = ke;
        pt.nemissions++;
    }

    // done
    return pt;
}

//
// Interpolate a value over an element
//
template <typename T>
static inline T interpolate_point(
    const T* vals, const int* points, int eid, const ym::vec3f& euv) {
    if (!vals) return T();
    auto& point = points[eid];
    return vals[point];
}

//
// Interpolate a value over an element
//
template <typename T>
static inline T interpolate_line(
    const T* vals, const ym::vec2i* lines, int eid, const ym::vec3f& euv) {
    if (!vals) return T();
    auto& line = lines[eid];
    return vals[line.x] * euv.x + vals[line.y] * euv.y;
}

//
// Interpolate a value over an element
//
template <typename T>
static inline T interpolate_triangle(
    const T* vals, const ym::vec3i* triangles, int eid, const ym::vec3f& euv) {
    if (!vals) return T();
    auto& triangle = triangles[eid];
    return vals[triangle.x] * euv.x + vals[triangle.y] * euv.y +
           vals[triangle.z] * euv.z;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static point eval_shapepoint(
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
    auto pos = ym::zero3f, norm = ym::zero3f;
    auto color = ym::zero4f;
    auto texcoord = ym::zero2f;
    auto tangsp = ym::zero4f;
    if (shp->points) {
        pos = interpolate_point(shp->pos, shp->points, eid, euv);
        norm =
            ym::normalize(interpolate_point(shp->norm, shp->points, eid, euv));
        texcoord = interpolate_point(shp->texcoord, shp->points, eid, euv);
        color = interpolate_point(shp->color, shp->points, eid, euv);
    } else if (shp->lines) {
        pos = interpolate_line(shp->pos, shp->lines, eid, euv);
        norm = ym::normalize(interpolate_line(shp->norm, shp->lines, eid, euv));
        texcoord = interpolate_line(shp->texcoord, shp->lines, eid, euv);
        color = interpolate_line(shp->color, shp->lines, eid, euv);
    } else if (shp->triangles) {
        pos = interpolate_triangle(shp->pos, shp->triangles, eid, euv);
        norm = interpolate_triangle(shp->norm, shp->triangles, eid, euv);
        texcoord =
            interpolate_triangle(shp->texcoord, shp->triangles, eid, euv);
        color = interpolate_triangle(shp->color, shp->triangles, eid, euv);
        tangsp = interpolate_triangle(shp->tangsp, shp->triangles, eid, euv);
    }

    // handle normal map
    if (shp->texcoord && shp->tangsp && shp->triangles && mat->norm_txt) {
        auto txt = eval_texture(mat->norm_txt, texcoord, false).xyz() * 2.0f -
                   ym::vec3f{1, 1, 1};
        auto ntxt = ym::normalize(ym::vec3f{txt.x, -txt.y, txt.z});
        auto frame = ym::make_frame3_fromzx(
            {0, 0, 0}, norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp[3];
        norm = ym::transform_direction(frame, ntxt);
    }

    // save texcoord
    pt.texcoord = texcoord;

    // creating frame
    pt.frame = ym::make_frame3_fromz(pos, norm);

    // transform to world space
    pt.frame.o = ym::transform_point(ist->frame, pt.frame.o);
    pt.frame.z = ym::transform_direction(ist->frame, pt.frame.z);
    pt.frame.x = ym::transform_direction(ist->frame, pt.frame.x);
    pt.frame.y = ym::transform_direction(ist->frame, pt.frame.y);

    // correct for doulbe sided
    if (mat->double_sided && ym::dot(pt.frame.z, pt.wo) < 0) {
        pt.frame.z = -pt.frame.z;
        pt.frame.x = -pt.frame.x;
    }

    // handle color
    auto kx_scale = ym::vec4f{1, 1, 1, 1};
    if (shp->color) kx_scale *= color;

    // handle occlusion
    if (shp->texcoord && mat->occ_txt)
        kx_scale.xyz() *= eval_texture(mat->occ_txt, texcoord).xyz();

    // sample emission
    auto ke = mat->ke * kx_scale.xyz();
    if (ke != ym::zero3f) {
        pt.emissions[0].type = emission_type::diffuse;
        pt.emissions[0].ke = ke;
        pt.nemissions++;
    }

    // sample reflectance
    switch (mat->rtype) {
        case reflectance_type::none: pt.op = 1; break;
        case reflectance_type::matte: {
            auto kd = ym::vec4f{mat->matte.kd, mat->matte.op} * kx_scale;
            if (shp->texcoord && mat->matte.kd_txt)
                kd *= eval_texture(mat->matte.kd_txt, texcoord);
            if (shp->texcoord && mat->matte.op_txt)
                kd.w *= eval_texture(mat->matte.op_txt, texcoord).x;
            pt.op = kd.w;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_lambert;
            pt.brdfs[pt.nbrdfs].rho = kd.xyz();
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
        } break;
        case reflectance_type::microfacet: {
            auto kd =
                ym::vec4f{mat->microfacet.kd, mat->microfacet.op} * kx_scale;
            auto ks = ym::vec4f{mat->microfacet.ks, mat->microfacet.rs} *
                      ym::vec4f{kx_scale.xyz(), 1};
            auto kt = ym::vec4f{mat->microfacet.kt, mat->microfacet.rs} *
                      ym::vec4f{kx_scale.xyz(), 1};
            if (shp->texcoord && mat->microfacet.kd_txt)
                kd *= eval_texture(mat->microfacet.kd_txt, texcoord);
            if (shp->texcoord && mat->microfacet.op_txt)
                kd.w *= eval_texture(mat->microfacet.op_txt, texcoord).x;
            if (shp->texcoord && mat->microfacet.ks_txt)
                ks.xyz() *=
                    eval_texture(mat->microfacet.ks_txt, texcoord).xyz();
            if (shp->texcoord && mat->microfacet.kt_txt)
                kt.xyz() *=
                    eval_texture(mat->microfacet.kt_txt, texcoord).xyz();
            pt.op = kd.w;
            pt.brdfs[pt.nbrdfs].type = brdf_type::refraction_ggx;
            pt.brdfs[pt.nbrdfs].rho = kt.xyz();
            pt.brdfs[pt.nbrdfs].roughness = kt.w;
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_lambert;
            pt.brdfs[pt.nbrdfs].rho = kd.xyz();
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_ggx;
            pt.brdfs[pt.nbrdfs].rho = ks.xyz();
            pt.brdfs[pt.nbrdfs].roughness = ks.w;
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
        } break;
        case reflectance_type::gltf_metallic_roughness: {
            auto kb =
                ym::vec4f{mat->metalrough.kb, mat->metalrough.op} * kx_scale;
            auto km = ym::vec2f{mat->metalrough.km, mat->metalrough.rs};
            if (shp->texcoord && mat->metalrough.kb_txt)
                kb *= eval_texture(mat->metalrough.kb_txt, texcoord);
            if (shp->texcoord && mat->metalrough.km_txt) {
                auto km_txt = eval_texture(mat->metalrough.km_txt, texcoord);
                km.x *= km_txt.y;
                km.y *= km_txt.z;
            }
            pt.op = kb.w;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_lambert;
            pt.brdfs[pt.nbrdfs].rho = kb.xyz() * (1 - km.x);
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_ggx;
            pt.brdfs[pt.nbrdfs].rho =
                kb.xyz() * km.x + ym::vec3f{0.04f, 0.04f, 0.04f} * (1 - km.x);
            pt.brdfs[pt.nbrdfs].roughness = km.y * km.y;
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f &&
                pt.brdfs[pt.nbrdfs].roughness < 0.999f)
                pt.nbrdfs++;
        } break;
        case reflectance_type::gltf_specular_glossiness: {
            auto kd =
                ym::vec4f{mat->specgloss.kd, mat->specgloss.op} * kx_scale;
            auto ks = ym::vec4f{mat->specgloss.ks, mat->specgloss.rs} *
                      ym::vec4f{kx_scale.xyz(), 1};
            if (shp->texcoord && mat->specgloss.kd_txt)
                kd *= eval_texture(mat->specgloss.kd_txt, texcoord);
            if (shp->texcoord && mat->specgloss.ks_txt)
                ks *= eval_texture(mat->specgloss.ks_txt, texcoord);
            pt.op = kd.w;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_lambert;
            pt.brdfs[pt.nbrdfs].rho = kd.xyz();
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_ggx;
            pt.brdfs[pt.nbrdfs].rho = ks.xyz();
            pt.brdfs[pt.nbrdfs].roughness = (1 - ks.w) * (1 - ks.w);
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
        } break;
        case reflectance_type::thin_glass: {
            auto ks = ym::vec4f{mat->thin_glass.ks, mat->microfacet.rs} *
                      ym::vec4f{kx_scale.xyz(), 1};
            auto kt = ym::vec4f{mat->thin_glass.kt, mat->microfacet.rs} *
                      ym::vec4f{kx_scale.xyz(), 1};
            if (shp->texcoord && mat->thin_glass.ks_txt)
                ks.xyz() *=
                    eval_texture(mat->thin_glass.ks_txt, texcoord).xyz();
            if (shp->texcoord && mat->thin_glass.kt_txt)
                kt.xyz() *=
                    eval_texture(mat->thin_glass.kt_txt, texcoord).xyz();
            pt.op = 1;
            pt.brdfs[pt.nbrdfs].type = brdf_type::transparent;
            pt.brdfs[pt.nbrdfs].rho = kt.xyz();
            pt.brdfs[pt.nbrdfs].roughness = 0;
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
            pt.brdfs[pt.nbrdfs].type = brdf_type::reflection_ggx;
            pt.brdfs[pt.nbrdfs].rho = ks.xyz();
            pt.brdfs[pt.nbrdfs].roughness = ks.w;
            if (pt.brdfs[pt.nbrdfs].rho != ym::zero3f) pt.nbrdfs++;
        } break;
    }

    // correct for different point types
    if (shp->points) {
        for (auto lid = 0; lid < pt.nbrdfs; lid++)
            pt.brdfs[lid].type = brdf_type::point_diffuse;
        for (auto lid = 0; lid < pt.nemissions; lid++)
            pt.emissions[lid].type = emission_type::point;
    } else if (shp->lines) {
        for (auto lid = 0; lid < pt.nbrdfs; lid++) {
            if (pt.brdfs[lid].type == brdf_type::reflection_lambert)
                pt.brdfs[lid].type = brdf_type::kajiya_kay_diff;
            else
                pt.brdfs[lid].type = brdf_type::kajiya_kay_spec;
        }
        for (auto lid = 0; lid < pt.nemissions; lid++)
            pt.emissions[lid].type = emission_type::line;
    }

    // correct for opacity
    if (pt.op < 1) {
        for (auto lid = 0; lid < pt.nemissions; lid++)
            pt.emissions[lid].ke *= pt.op;
        for (auto lid = 0; lid < pt.nbrdfs; lid++) pt.brdfs[lid].rho *= pt.op;
        pt.brdfs[pt.nbrdfs].type = brdf_type::transparent;
        pt.brdfs[pt.nbrdfs].rho = {1 - pt.op, 1 - pt.op, 1 - pt.op};
        pt.nbrdfs += 1;
    }

    // set whole albedo
    pt.rho = ym::zero3f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) pt.rho += pt.brdfs[lid].rho;

    return pt;
}

//
// Sample weight for a light point.
//
static float weight_light(const point& lpt, const point& pt) {
    if (lpt.no_emission()) return 0;
    // support only one lobe for now
    auto emission = lpt.emissions[0];
    switch (emission.type) {
        case emission_type::env: {
            return 4 * ym::pif;
        } break;
        case emission_type::point: {
            auto d = ym::dist(lpt.frame.o, pt.frame.o);
            return lpt.ist->shp->area / (d * d);
        } break;
        case emission_type::line: {
            assert(false);
            return 0;
        } break;
        case emission_type::diffuse: {
            auto d = ym::dist(lpt.frame.o, pt.frame.o);
            return lpt.ist->shp->area * std::abs(ym::dot(lpt.frame.z, lpt.wo)) /
                   (d * d);
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
static point sample_light(
    const light* lgt, const point& pt, float rne, const ym::vec2f& rn) {
    if (lgt->ist) {
        auto shp = lgt->ist->shp;
        auto eid = 0;
        auto euv = ym::zero3f;
        if (shp->triangles) {
            std::tie(eid, euv) =
                ym::sample_triangles(shp->nelems, shp->cdf.data(), rne, rn);
        } else if (shp->lines) {
            std::tie(eid, (ym::vec2f&)euv) =
                ym::sample_lines(shp->nelems, shp->cdf.data(), rne, rn.x);
        } else if (shp->points) {
            eid = ym::clamp(0, shp->nelems - 1, (int)(rne * shp->nelems));
            euv = {1, 0, 0};
        } else {
            assert(false);
        }
        auto lpt = eval_shapepoint(lgt->ist, eid, euv, ym::zero3f);
        lpt.wo = ym::normalize(pt.frame.o - lpt.frame.o);
        return lpt;
    } else if (lgt->env) {
        auto z = -1 + 2 * rn.y;
        auto rr = std::sqrt(ym::clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * ym::pif * rn.x;
        auto wo = ym::vec3f{std::cos(phi) * rr, z, std::sin(phi) * rr};
        auto lpt = eval_envpoint(lgt->env, wo);
        return lpt;
    } else {
        assert(false);
        return {};
    }
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f offset_ray(
    const point& pt, const ym::vec3f& w, const trace_params& params) {
    if (dot(w, pt.frame.z) > 0) {
        return ym::ray3f(
            pt.frame.o + pt.frame.z * params.ray_eps, w, params.ray_eps);
    } else {
        return ym::ray3f(
            pt.frame.o - pt.frame.z * params.ray_eps, w, params.ray_eps);
    }
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f offset_ray(
    const point& pt, const point& pt2, const trace_params& params) {
    auto ray_dist = (!pt2.env) ? ym::dist(pt.frame.o, pt2.frame.o) : FLT_MAX;
    if (dot(pt2.frame.o - pt.frame.o, pt.frame.z) > 0) {
        return ym::ray3f(pt.frame.o + pt.frame.z * params.ray_eps, -pt2.wo,
            params.ray_eps, ray_dist - 2 * params.ray_eps);
    } else {
        return ym::ray3f(pt.frame.o - pt.frame.z * params.ray_eps, -pt2.wo,
            params.ray_eps, ray_dist - 2 * params.ray_eps);
    }
}

//
// Intersects a ray with the scn and return the point (or env point).
//
static point intersect_scene(const scene* scn, const ym::ray3f& ray) {
    auto isec = scn->intersect_first(ray);
    if (isec) {
        return eval_shapepoint(
            scn->instances[isec.iid], isec.eid, isec.euv, -ray.d);
    } else if (!scn->environments.empty()) {
        return eval_envpoint(scn->environments[0], -ray.d);
    } else {
        return {};
    }
}

//
// Transparecy
//
static ym::vec3f eval_transparency(const point& pt) {
    auto kt = ym::zero3f;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        auto brdf = pt.brdfs[lid];
        switch (brdf.type) {
            case brdf_type::transparent:
            case brdf_type::refraction_ggx: kt += brdf.rho; break;
            default: break;
        }
    }
    return kt;
}

//
// Test occlusion
//
static ym::vec3f eval_transmission(const scene* scn, const point& pt,
    const point& lpt, const trace_params& params) {
    if (scn->shadow_transmission) {
        auto cpt = pt;
        auto weight = ym::vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            cpt = intersect_scene(scn, offset_ray(cpt, lpt, params));
            if (!cpt.ist) break;
            weight *= eval_transparency(cpt);
            if (weight == ym::zero3f) break;
        }
        return weight;
    } else {
        auto shadow_ray = offset_ray(pt, lpt, params);
        return (scn->intersect_any(shadow_ray)) ? ym::zero3f :
                                                  ym::vec3f{1, 1, 1};
    }
}

//
// Mis weight
//
static inline float weight_mis(float w0, float w1) {
    if (!w0 || !w1) return 1;
    return (1 / w0) / (1 / w0 + 1 / w1);
}

//
// Recursive path tracing.
//
static ym::vec3f shade_pathtrace(const scene* scn, const point& pt_,
    sampler* smp, const trace_params& params) {
    // make a copy
    auto pt = pt_;

    // emission
    auto l = eval_emission(pt);
    if (pt.no_reflectance() || scn->lights.empty()) return l;

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt);

        // direct – light
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto lw = weight_light(lpt, pt) * (float)scn->lights.size();
        auto lke = eval_emission(lpt);
        auto lbc = eval_brdfcos(pt, -lpt.wo);
        auto lld = lke * lbc * lw;
        if (lld != ym::zero3f) {
            l += weight * lld * eval_transmission(scn, pt, lpt, params) *
                 weight_mis(lw, weight_brdfcos(pt, -lpt.wo));
        }

        // direct – brdf
        auto bpt = intersect_scene(
            scn, offset_ray(pt,
                     sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp)),
                     params));
        auto bw = weight_brdfcos(pt, -bpt.wo);
        auto bke = eval_emission(bpt);
        auto bbc = eval_brdfcos(pt, -bpt.wo);
        auto bld = bke * bbc * bw;
        if (bld != ym::zero3f) {
            l += weight * bld * weight_mis(bw, weight_light(bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;
        if (bpt.no_reflectance()) break;

        // continue path
        weight *= eval_brdfcos(pt, -bpt.wo) * weight_brdfcos(pt, -bpt.wo);
        if (weight == ym::zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - std::min(std::max(std::max(pt.rho.x, pt.rho.y),
                                              pt.rho.z),
                                     0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        emission = false;
    }

    return l;
}

//
// Recursive path tracing.
//
static ym::vec3f shade_pathtrace_std(const scene* scn, const point& pt_,
    sampler* smp, const trace_params& params) {
    // amke a copy
    auto pt = pt_;

    // emission
    auto l = eval_emission(pt);
    if (pt.no_reflectance()) return l;

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt);

        // direct
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != ym::zero3f) {
            auto shadow_ray = offset_ray(pt, lpt, params);
            if (!scn->intersect_any(shadow_ray)) l += weight * ld;
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - std::min(std::max(std::max(pt.rho.x, pt.rho.y),
                                              pt.rho.z),
                                     0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        {
            auto wi =
                sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp));
            weight *= eval_brdfcos(pt, wi) * weight_brdfcos(pt, wi);
            if (weight == ym::zero3f) break;

            pt = intersect_scene(scn, offset_ray(pt, wi, params));
            emission = false;
            if (pt.no_reflectance()) break;
        }
    }

    return l;
}

//
// Recursive path tracing.
//
static ym::vec3f shade_pathtrace_hack(const scene* scn, const point& pt_,
    sampler* smp, const trace_params& params) {
    // make a copy
    auto pt = pt_;

    // emission
    auto l = eval_emission(pt);
    if (pt.no_reflectance()) return l;

    // trace path
    auto weight = ym::vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // direct
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != ym::zero3f) {
            auto shadow_ray = offset_ray(pt, lpt, params);
            if (!scn->intersect_any(shadow_ray)) l += weight * ld;
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - std::min(std::max(std::max(pt.rho.x, pt.rho.y),
                                              pt.rho.z),
                                     0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        {
            auto wi =
                sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp));
            weight *= eval_brdfcos(pt, wi) * weight_brdfcos(pt, wi);
            if (weight == ym::zero3f) break;

            pt = intersect_scene(scn, offset_ray(pt, wi, params));
            if (pt.no_reflectance()) break;
        }
    }

    return l;
}

//
// Direct illumination.
//
static ym::vec3f shade_direct(const scene* scn, const point& pt, int bounce,
    sampler* smp, const trace_params& params) {
    // emission
    auto l = eval_emission(pt);
    if (pt.no_reflectance()) return l;

    // ambient
    l += params.amb * pt.rho;

    // direct
    for (auto& lgt : scn->lights) {
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt);
        if (ld == ym::zero3f) continue;
        l += ld * eval_transmission(scn, pt, lpt, params);
    }

    // exit if needed
    if (bounce >= params.max_depth) return l;

    // opacity
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        auto& brdf = pt.brdfs[lid];
        if (brdf.type == brdf_type::transparent) {
            auto ray = offset_ray(pt, -pt.wo, params);
            l += brdf.rho * shade_direct(scn, intersect_scene(scn, ray),
                                bounce + 1, smp, params);
        }
    }

    // done
    return l;
}

//
// Direct illumination.
//
static ym::vec3f shade_direct(const scene* scn, const point& pt, sampler* smp,
    const trace_params& params) {
    return shade_direct(scn, pt, 0, smp, params);
}

//
// Eyelight for quick previewing.
//
static ym::vec3f shade_eyelight(const scene* scn, const point& pt, int bounce,
    sampler* smp, const trace_params& params) {
    // emission
    auto l = eval_emission(pt);
    if (pt.no_reflectance()) return l;

    // brdf*light
    l += eval_brdfcos(pt, pt.wo) * ym::pif;

    // opacity
    if (bounce >= params.max_depth) return l;
    for (auto lid = 0; lid < pt.nbrdfs; lid++) {
        auto& brdf = pt.brdfs[lid];
        if (brdf.type == brdf_type::transparent) {
            auto ray = offset_ray(pt, -pt.wo, params);
            l += brdf.rho * shade_eyelight(scn, intersect_scene(scn, ray),
                                bounce + 1, smp, params);
        }
    }

    // done
    return l;
}

//
// Eyelight for quick previewing.
//
static ym::vec3f shade_eyelight(const scene* scn, const point& pt, sampler* smp,
    const trace_params& params) {
    return shade_eyelight(scn, pt, 0, smp, params);
}

//
// Debug previewing.
//
static ym::vec3f shade_debug_normal(const scene* scn, const point& pt,
    sampler* smp, const trace_params& params) {
    return pt.frame.z * 0.5f + ym::vec3f{0.5f, 0.5f, 0.5f};
}

//
// Debug previewing.
//
static ym::vec3f shade_debug_albedo(const scene* scn, const point& pt,
    sampler* smp, const trace_params& params) {
    return pt.rho;
}

//
// Debug previewing.
//
static ym::vec3f shade_debug_texcoord(const scene* scn, const point& pt,
    sampler* smp, const trace_params& params) {
    return {pt.texcoord.x, pt.texcoord.y, 0};
}

//
// Shader function callback.
//
using shade_fn = ym::vec3f (*)(const scene* scn, const point& pt, sampler* smp,
    const trace_params& params);

//
// Renders a block of pixels. Public API, see above.
//
void trace_block(const scene* scn, ym::vec4f* img, int block_x, int block_y,
    int block_width, int block_height, int samples_min, int samples_max,
    const trace_params& params) {
    auto cam = scn->cameras[params.camera_id];
    shade_fn shade;
    switch (params.stype) {
        case shader_type::eyelight: shade = shade_eyelight; break;
        case shader_type::direct: shade = shade_direct; break;
        case shader_type::pathtrace: shade = shade_pathtrace; break;
        default: assert(false); return;
    }
    for (auto j = block_y; j < block_y + block_height; j++) {
        for (auto i = block_x; i < block_x + block_width; i++) {
            auto lp = ym::zero4f;
            for (auto s = samples_min; s < samples_max; s++) {
                auto smp = make_sampler(i, j, s, params.nsamples, params.rtype);
                auto rn = sample_next2f(&smp);
                auto uv = ym::vec2f{
                    (i + rn.x) / params.width, 1 - (j + rn.y) / params.height};
                auto ray = eval_camera(cam, uv, sample_next2f(&smp));
                auto pt = intersect_scene(scn, ray);
                if (!pt.ist || params.envmap_invisible) continue;
                auto l = shade(scn, pt, &smp, params);
                if (!ym::isfinite(l)) {
                    if (scn->log_error) scn->log_error("NaN detected");
                    continue;
                }
                if (params.pixel_clamp > 0)
                    l = ym::clamplen(l, params.pixel_clamp);
                lp += {l, 1};
            }
            img[j * params.width + i] = lp / (float)(samples_max - samples_min);
        }
    }
}

// triangle filter (public domain from stb_image_resize)
inline float filter_triangle(float x) {
    x = (float)fabs(x);

    if (x <= 1.0f)
        return 1 - x;
    else
        return 0;
}

// cubic filter (public domain from stb_image_resize)
inline float filter_cubic(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return (4 + x * x * (3 * x - 6)) / 6;
    else if (x < 2.0f)
        return (8 + x * (-12 + x * (6 - x))) / 6;
    else
        return 0.0f;
}

// catmull-rom filter (public domain from stb_image_resize)
inline float filter_catmullrom(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return 1 - x * x * (2.5f - 1.5f * x);
    else if (x < 2.0f)
        return 2 - x * (4 + x * (0.5f * x - 2.5f));
    else
        return 0.0f;
}

// mitchell filter (public domain from stb_image_resize)
inline float filter_mitchell(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return (16 + x * x * (21 * x - 36)) / 18;
    else if (x < 2.0f)
        return (32 + x * (-60 + x * (36 - 7 * x))) / 18;
    else
        return 0.0f;
}

// filter function
using filter_fn = float (*)(float);

//
// state for progressive rendering and denoising
//
struct trace_state {
    // rendered image
    ym::image4f img;

    // progressive rendering buffers
    ym::image4f acc;
    ym::imagef weight;

    // auxiliary buffers
    ym::image4f norm, albedo, depth;

    // progressive state
    int cur_sample = 0;
    std::vector<ym::bbox2i> blocks;

    // pool
    yu::concurrent::thread_pool* pool = nullptr;
    /// lock for access to image
    std::mutex image_mutex;

    // render scene
    const scene* scn = nullptr;
    // render options
    trace_params params = {};

    // render function
    shade_fn shade = nullptr;
    // render camera
    const camera* cam = nullptr;
    // filter function
    filter_fn filter = nullptr;
    // filter size
    int filter_size = 0;

    // cleanup
    ~trace_state() {
        if (pool) {
            yu::concurrent::clear_pool(pool);
            yu::concurrent::free_pool(pool);
        }
    }
};

//
// Make image blocks
//
std::vector<ym::bbox2i> make_blocks(int w, int h, int bs) {
    std::vector<ym::bbox2i> blocks;
    for (int j = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs) {
            blocks.push_back(
                {{i, j}, {ym::min(i + bs, w), ym::min(j + bs, h)}});
        }
    }
    return blocks;
}

//
// Initialize state
//
trace_state* make_state() { return new trace_state(); }

//
// Initialize state
//
void init_state(
    trace_state* state, const scene* scn, const trace_params& params) {
    if (state->pool) {
        if (params.parallel)
            yu::concurrent::clear_pool(state->pool);
        else
            yu::concurrent::free_pool(state->pool);
    } else {
        if (params.parallel) state->pool = yu::concurrent::make_pool();
    }
    state->img = ym::image4f(params.width, params.height);
    state->acc = ym::image4f(params.width, params.height);
    state->weight = ym::imagef(params.width, params.height);
    if (params.aux_buffers) {
        state->norm = ym::image4f(params.width, params.height);
        state->albedo = ym::image4f(params.width, params.height);
        state->depth = ym::image4f(params.width, params.height);
    } else {
        state->norm = ym::image4f();
        state->albedo = ym::image4f();
        state->depth = ym::image4f();
    }
    state->cur_sample = 0;
    state->blocks = make_blocks(params.width, params.height, 32);
    state->scn = scn;
    state->params = params;

    state->cam = scn->cameras[params.camera_id];

    switch (params.stype) {
        case shader_type::eyelight: state->shade = shade_eyelight; break;
        case shader_type::direct: state->shade = shade_direct; break;
        case shader_type::pathtrace: state->shade = shade_pathtrace; break;
        case shader_type::debug_normal:
            state->shade = shade_debug_normal;
            break;
        case shader_type::debug_albedo:
            state->shade = shade_debug_albedo;
            break;
        case shader_type::debug_texcoord:
            state->shade = shade_debug_texcoord;
            break;
        default: assert(false); return;
    }

    switch (params.ftype) {
        case filter_type::box:
            state->filter = nullptr;
            state->filter_size = 0;
            break;
        case filter_type::triangle:
            state->filter = filter_triangle;
            state->filter_size = 1;
            break;
        case filter_type::cubic:
            state->filter = filter_cubic;
            state->filter_size = 2;
            break;
        case filter_type::catmull_rom:
            state->filter = filter_catmullrom;
            state->filter_size = 2;
            break;
        case filter_type::mitchell:
            state->filter = filter_mitchell;
            state->filter_size = 2;
            break;
        default: assert(false); return;
    }
}

//
// Grabs the image from the state
//
ym::image4f& get_traced_image(trace_state* state) { return state->img; }

//
// Grabs the image from the state
//
void get_aux_buffers(const trace_state* state, ym::image4f& norm,
    ym::image4f& albedo, ym::image4f& depth) {
    if (!state->params.aux_buffers) return;
    norm.resize(state->norm.width(), state->norm.height());
    albedo.resize(state->albedo.width(), state->albedo.height());
    depth.resize(state->depth.width(), state->depth.height());
    for (auto j = 0; j < norm.height(); j++) {
        for (auto i = 0; i < norm.width(); i++) {
            norm[{i, j}] = state->norm[{i, j}] / state->weight[{i, j}];
            norm[{i, j}].xyz() =
                normalize(norm[{i, j}].xyz()) * 0.5f + ym::vec3f{0.5, 0.5, 0.5};
            albedo[{i, j}] = state->albedo[{i, j}] / state->weight[{i, j}];
            depth[{i, j}] = state->depth[{i, j}] / (10 * state->weight[{i, j}]);
        }
    }
}

//
// Grt the current sample count
//
int get_cur_sample(const trace_state* state) { return state->cur_sample; }

//
// Trace a single sample
//
void trace_sample(ytrace::trace_state* state, int i, int j, int s, ym::vec3f& l,
    ytrace::point& pt, ym::vec2f& rn) {
    auto& params = state->params;
    auto smp = make_sampler(i, j, s, params.nsamples, params.rtype);
    rn = sample_next2f(&smp);
    auto uv =
        ym::vec2f{(i + rn.x) / params.width, 1 - (j + rn.y) / params.height};
    auto ray = eval_camera(state->cam, uv, sample_next2f(&smp));
    pt = intersect_scene(state->scn, ray);
    if (!pt.ist || params.envmap_invisible) return;
    l = state->shade(state->scn, pt, &smp, state->params);
    if (!ym::isfinite(l)) {
        if (state->scn->log_error) state->scn->log_error("NaN detected");
        return;
    }
    if (params.pixel_clamp > 0) l = ym::clamplen(l, params.pixel_clamp);
}

//
// Trace a block of samples
//
void trace_block_box(
    trace_state* state, int block_idx, int samples_min, int samples_max) {
    auto& block = state->blocks[block_idx];
    for (auto j = block.min.y; j < block.max.y; j++) {
        for (auto i = block.min.x; i < block.max.x; i++) {
            for (auto s = samples_min; s < samples_max; s++) {
                ytrace::point pt;
                auto l = ym::zero3f;
                auto uv = ym::zero2f;
                trace_sample(state, i, j, s, l, pt, uv);
                state->acc[{i, j}] += {l, 1};
                state->weight[{i, j}] += 1;
                state->img[{i, j}] = state->acc[{i, j}] / state->weight[{i, j}];
                if (state->params.aux_buffers && pt.ist) {
                    state->norm[{i, j}] += {pt.frame.z, 1};
                    state->albedo[{i, j}] += {pt.rho, 1};
                    auto d = length(pt.frame.o - state->cam->frame.o);
                    state->depth[{i, j}] += {d, d, d, 1};
                }
            }
        }
    }
}

//
// Trace a block of samples
//
void trace_block_filtered(
    trace_state* state, int block_idx, int samples_min, int samples_max) {
    static constexpr const int pad = 2;
    auto& block = state->blocks[block_idx];
    auto block_size = ym::diagonal(block);
    auto acc_buffer =
        ym::image4f(block_size.x + pad * 2, block_size.y + pad * 2);
    auto weight_buffer =
        ym::imagef(block_size.x + pad * 2, block_size.y + pad * 2);
    for (auto j = block.min.y; j < block.max.y; j++) {
        for (auto i = block.min.x; i < block.max.x; i++) {
            for (auto s = samples_min; s < samples_max; s++) {
                ytrace::point pt;
                auto l = ym::zero3f;
                auto uv = ym::zero2f;
                trace_sample(state, i, j, s, l, pt, uv);
                if (state->filter) {
                    auto bi = i - block.min.x, bj = j - block.min.y;
                    for (auto fj = -state->filter_size;
                         fj <= state->filter_size; fj++) {
                        for (auto fi = -state->filter_size;
                             fi <= state->filter_size; fi++) {
                            auto w = filter_triangle(fi - uv.x + 0.5f) *
                                     filter_triangle(fj - uv.y + 0.5f);
                            acc_buffer[{bi + fi + pad, bj + fj + pad}] +=
                                {l * w, w};
                            weight_buffer[{bi + fi + pad, bj + fj + pad}] += w;
                        }
                    }
                } else {
                    auto bi = i - block.min.x, bj = j - block.min.y;
                    acc_buffer[{bi + pad, bj + pad}] += {l, 1};
                    weight_buffer[{bi + pad, bj + pad}] += 1;
                }
            }
        }
    }
    if (state->filter) {
        std::unique_lock<std::mutex> lock_guard(state->image_mutex);
        auto width = state->acc.width(), height = state->acc.height();
        for (auto j = ym::max(block.min.y - state->filter_size, 0);
             j < ym::min(block.max.y + state->filter_size, height); j++) {
            for (auto i = ym::max(block.min.x - state->filter_size, 0);
                 i < ym::min(block.max.x + state->filter_size, width); i++) {
                auto bi = i - block.min.x, bj = j - block.min.y;
                state->acc[{i, j}] += acc_buffer[{bi + pad, bj + pad}];
                state->weight[{i, j}] += weight_buffer[{bi + pad, bj + pad}];
                state->img[{i, j}] = state->acc[{i, j}] / state->weight[{i, j}];
            }
        }
    } else {
        for (auto j = block.min.y; j < block.max.y; j++) {
            for (auto i = block.min.x; i < block.max.x; i++) {
                auto bi = i - block.min.x, bj = j - block.min.y;
                state->acc[{i, j}] += acc_buffer[{bi + pad, bj + pad}];
                state->weight[{i, j}] += weight_buffer[{bi + pad, bj + pad}];
                state->img[{i, j}] = state->acc[{i, j}] / state->weight[{i, j}];
            }
        }
    }
}

//
// Trace a block of samples
//
void trace_block(
    trace_state* state, int block_idx, int samples_min, int samples_max) {
    if (state->filter)
        return trace_block_filtered(state, block_idx, samples_min, samples_max);
    else
        return trace_block_box(state, block_idx, samples_min, samples_max);
}
//
// Clear state
//
void free_state(trace_state*& state) {
    if (state) delete state;
    state = nullptr;
}

//
// Trace a batch of samples.
//
bool trace_next_samples(trace_state* state, int nsamples) {
    if (state->cur_sample >= state->params.nsamples) return false;
    nsamples = ym::min(nsamples, state->params.nsamples - state->cur_sample);
    if (state->pool) {
        yu::concurrent::parallel_for(
            state->pool, (int)state->blocks.size(), [state, nsamples](int idx) {
                ytrace::trace_block(state, idx, state->cur_sample,
                    state->cur_sample + nsamples);
            });
    } else {
        for (auto idx = 0; idx < (int)state->blocks.size(); idx++) {
            ytrace::trace_block(
                state, idx, state->cur_sample, state->cur_sample + nsamples);
        }
    }
    state->cur_sample += nsamples;
    return true;
}

//
// Starts an anyncrhounous renderer with a maximum of 256 samples.
//
void trace_async_start(trace_state* state) {
    for (auto sample = 0; sample < state->params.nsamples; sample++) {
        for (auto block_idx = 0; block_idx < state->blocks.size();
             block_idx++) {
            yu::concurrent::run_async(
                state->pool, [state, sample, block_idx]() {
                    ytrace::trace_block(state, block_idx, sample, sample + 1);
                    if (block_idx == state->blocks.size() - 1) {
                        state->cur_sample++;
                    }
                });
        }
    }
}

//
// Stop the asynchronous renderer.
//
void trace_async_stop(trace_state* state) {
    if (!state->pool) return;
    yu::concurrent::clear_pool(state->pool);
}

}  // namespace ytrace

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif
