//
// Implementation for Yocto/Scene.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "yocto_scene.h"
#include "yocto_bvh.h"
#include "yocto_image.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include <unordered_map>

#if YGL_OBJ
#include "yocto_obj.h"
#endif

#if YGL_GLTF
#include "yocto_gltf.h"
#endif

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {
shape::~shape() {
    if (bvh) delete bvh;
}
scene::~scene() {
    for (auto v : shapes) delete v;
    for (auto v : instances) delete v;
    for (auto v : materials) delete v;
    for (auto v : textures) delete v;
    for (auto v : cameras) delete v;
    for (auto v : environments) delete v;
    for (auto v : nodes) delete v;
    for (auto v : animations) delete v;
    if (bvh) delete bvh;
}
}  // namespace ygl

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute shape normals. Supports only non-facevarying shapes.
void update_normals(shape* shp) {
    if (!shp->triangles.empty()) {
        compute_normals(shp->triangles, shp->pos, shp->norm);
    } else if (!shp->lines.empty()) {
        compute_tangents(shp->lines, shp->pos, shp->norm);
    } else {
        shp->norm.assign(shp->pos.size(), {0, 0, 1});
    }
}

// Computes a shape bounding box.
void update_bbox(shape* shp) {
    shp->bbox = invalid_bbox3f;
    for (auto p : shp->pos) shp->bbox += p;
}

// Updates the scene and scene's instances bounding boxes
void update_bbox(scene* scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bbox(shp);
    }
    scn->bbox = invalid_bbox3f;
    for (auto ist : scn->instances) {
        ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);
        scn->bbox += ist->bbox;
    }
}

// Update animation transforms
void update_transforms(
    const animation* anm, float time, const std::string& anim_group) {
    if (anim_group != "" && anim_group != anm->group) return;

    if (!anm->translation.empty()) {
        auto val = vec3f{0, 0, 0};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->translation, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->translation, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->translation, time);
                break;
            default: throw std::runtime_error("should not have been here");
        }
        for (auto target : anm->targets) target->translation = val;
    }
    if (!anm->rotation.empty()) {
        auto val = vec4f{0, 0, 0, 1};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->rotation, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->rotation, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->rotation, time);
                break;
        }
        for (auto target : anm->targets) target->rotation = val;
    }
    if (!anm->scale.empty()) {
        auto val = vec3f{1, 1, 1};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->scale, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->scale, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->scale, time);
                break;
        }
        for (auto target : anm->targets) target->scale = val;
    }
}

// Update node transforms
void update_transforms(node* nde, const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame(nde->translation) *
                 rotation_frame(nde->rotation) * scaling_frame(nde->scale);
    if (nde->ist) nde->ist->frame = frame;
    if (nde->cam) nde->cam->frame = frame;
    if (nde->env) nde->env->frame = frame;
    for (auto child : nde->children) update_transforms(child, frame);
}

// Update node transforms
void update_transforms(scene* scn, float time, const std::string& anim_group) {
    for (auto agr : scn->animations) update_transforms(agr, time, anim_group);
    for (auto nde : scn->nodes) nde->children.clear();
    for (auto nde : scn->nodes)
        if (nde->parent) nde->parent->children.push_back(nde);
    for (auto nde : scn->nodes)
        if (!nde->parent) update_transforms(nde);
}

// Compute animation range
vec2f compute_animation_range(const scene* scn, const std::string& anim_group) {
    if (scn->animations.empty()) return zero2f;
    auto range = vec2f{+flt_max, -flt_max};
    for (auto anm : scn->animations) {
        if (anim_group != "" && anm->group != anim_group) continue;
        range.x = min(range.x, anm->times.front());
        range.y = max(range.y, anm->times.back());
    }
    if (range.y < range.x) return zero2f;
    return range;
}

// Update lights.
void update_lights(scene* scn, bool do_shapes, bool do_environments) {
    if (do_shapes) {
        for (auto shp : scn->shapes) shp->elem_cdf.clear();
    }
    if (do_environments) {
        for (auto env : scn->environments) env->elem_cdf.clear();
    }
    scn->lights.clear();

    for (auto ist : scn->instances) {
        if (!ist->mat || ist->mat->ke == zero3f) continue;
        scn->lights.push_back(ist);
        if (ist->shp->elem_cdf.empty()) update_shape_cdf(ist->shp);
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        if (env->elem_cdf.empty()) update_environment_cdf(env);
    }
}

// Generate a distribution for sampling a shape uniformly based
// on area/length.
void update_shape_cdf(shape* shp) {
    shp->elem_cdf.clear();
    if (!shp->triangles.empty()) {
        shp->elem_cdf = sample_triangles_cdf(shp->triangles, shp->pos);
    } else if (!shp->lines.empty()) {
        shp->elem_cdf = sample_lines_cdf(shp->lines, shp->pos);
    } else if (!shp->pos.empty()) {
        shp->elem_cdf = sample_points_cdf(shp->pos.size());
    } else {
        throw std::runtime_error("empty shape not supported");
    }
}

// Update environment CDF for sampling.
void update_environment_cdf(environment* env) {
    env->elem_cdf.clear();
    auto txt = env->ke_txt.txt;
    if (!txt) return;
    env->elem_cdf.resize(txt->width * txt->height);
    if (!txt->ldr.empty()) {
        for (auto i = 0; i < env->elem_cdf.size(); i++) {
            auto th = (i / txt->width + 0.5f) * pi / txt->height;
            auto rgba = srgb_to_linear(txt->ldr[i]);
            env->elem_cdf[i] = max(rgba.x, max(rgba.y, rgba.z)) * sin(th);
            if (i) env->elem_cdf[i] += env->elem_cdf[i - 1];
        }
    } else if (!txt->hdr.empty()) {
        for (auto i = 0; i < env->elem_cdf.size(); i++) {
            auto th = (i / txt->width + 0.5f) * pi / txt->height;
            auto rgba = txt->hdr[i];
            env->elem_cdf[i] = max(rgba.x, max(rgba.y, rgba.z)) * sin(th);
            if (i) env->elem_cdf[i] += env->elem_cdf[i - 1];
        }
    } else {
        throw std::runtime_error("empty texture");
    }
}

// Build a shape BVH
void update_bvh(shape* shp, bool equalsize) {
    if (!shp->bvh) shp->bvh = new bvh_tree();
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    shp->bvh->lines = shp->lines;
    shp->bvh->triangles = shp->triangles;
    build_bvh(shp->bvh, equalsize);
}

// Build a scene BVH
void update_bvh(scene* scn, bool do_shapes, bool equalsize) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bvh(shp, equalsize);
    }

    // tree bvh
    if (!scn->bvh) scn->bvh = new bvh_tree();
    scn->bvh->ist_frames.resize(scn->instances.size());
    scn->bvh->ist_inv_frames.resize(scn->instances.size());
    scn->bvh->ist_bvhs.resize(scn->instances.size());
    for (auto i = 0; i < scn->instances.size(); i++) {
        auto ist = scn->instances[i];
        scn->bvh->ist_frames[i] = ist->frame;
        scn->bvh->ist_inv_frames[i] = inverse(ist->frame);
        scn->bvh->ist_bvhs[i] = ist->shp->bvh;
    }
    build_bvh(scn->bvh, equalsize);
}

// Refits a scene BVH
void refit_bvh(shape* shp) {
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    refit_bvh(shp->bvh);
}

// Refits a scene BVH
void refit_bvh(scene* scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) refit_bvh(shp);
    }
    scn->bvh->ist_frames.resize(scn->instances.size());
    scn->bvh->ist_inv_frames.resize(scn->instances.size());
    scn->bvh->ist_bvhs.resize(scn->instances.size());
    for (auto i = 0; i < scn->instances.size(); i++) {
        auto ist = scn->instances[i];
        scn->bvh->ist_frames[i] = ist->frame;
        scn->bvh->ist_inv_frames[i] = inverse(ist->frame);
        scn->bvh->ist_bvhs[i] = ist->shp->bvh;
    }
    refit_bvh(scn->bvh);
}

// Add missing names and resolve duplicated names.
void add_missing_names(scene* scn) {
    auto fix_names = [](auto& vals, const std::string& base) {
        auto nmap = std::map<std::string, int>();
        for (auto val : vals) {
            if (val->name == "") val->name = base;
            if (nmap.find(val->name) == nmap.end()) {
                nmap[val->name] = 0;
            } else {
                nmap[val->name] += 1;
                val->name = val->name + "_" + std::to_string(nmap[val->name]);
            }
        }
    };
    fix_names(scn->cameras, "cam");
    fix_names(scn->shapes, "shp");
    fix_names(scn->textures, "txt");
    fix_names(scn->materials, "mat");
    fix_names(scn->environments, "env");
    fix_names(scn->nodes, "nde");
    fix_names(scn->animations, "anm");
}

// Add missing normals.
void add_missing_normals(scene* scn) {
    for (auto shp : scn->shapes) {
        if (shp->norm.empty()) update_normals(shp);
    }
}

// Add missing tangent space if needed.
void add_missing_tangent_space(scene* scn) {
    for (auto ist : scn->instances) {
        if (!ist->shp->tangsp.empty() || ist->shp->texcoord.empty()) continue;
        if (!ist->mat || (!ist->mat->norm_txt.txt && !ist->mat->bump_txt.txt))
            continue;
        if (!ist->shp->triangles.empty()) {
            if (ist->shp->norm.empty()) update_normals(ist->shp);
            compute_tangent_frames(ist->shp->triangles, ist->shp->pos,
                ist->shp->norm, ist->shp->texcoord, ist->shp->tangsp);
        } else {
            throw std::runtime_error("type not supported");
        }
    }
}

// add missing camera
void add_missing_camera(scene* scn) {
    if (!scn->cameras.empty()) return;
    update_bbox(scn, true);
    auto bbox = scn->bbox;
    auto bbox_center = (bbox.max + bbox.min) / 2.0f;
    auto bbox_size = bbox.max - bbox.min;
    auto bbox_msize = max(bbox_size.x, max(bbox_size.y, bbox_size.z));
    auto cam = new camera();
    cam->name = "<view>";
    auto camera_dir = vec3f{1, 0.4f, 1};
    auto from = camera_dir * bbox_msize + bbox_center;
    auto to = bbox_center;
    auto up = vec3f{0, 1, 0};
    cam->frame = lookat_frame(from, to, up);
    cam->ortho = false;
    cam->aspect = 16.0f / 9.0f;
    cam->yfov = 2 * atanf(0.5f);
    cam->aperture = 0;
    cam->focus = length(to - from);
    scn->cameras.push_back(cam);
}

// Checks for validity of the scene.
std::vector<std::string> validate(const scene* scn, bool skip_textures) {
    auto errs = std::vector<std::string>();
    auto check_names = [&errs](const auto& vals, const std::string& base) {
        auto used = std::map<std::string, int>();
        for (auto val : vals) used[val->name] += 1;
        for (auto& kv : used) {
            if (kv.first == "")
                errs.push_back("empty " + base + " name");
            else if (kv.second > 1)
                errs.push_back("duplicated " + base + " name " + kv.first);
        }
    };
    auto check_empty_textures = [&errs](const std::vector<texture*>& vals) {
        for (auto val : vals) {
            if (val->ldr.empty() && val->hdr.empty())
                errs.push_back("empty texture " + val->name);
        }
    };

    check_names(scn->cameras, "camera");
    check_names(scn->shapes, "shape");
    check_names(scn->textures, "texture");
    check_names(scn->materials, "material");
    check_names(scn->environments, "environment");
    check_names(scn->nodes, "node");
    check_names(scn->animations, "animation");
    if (!skip_textures) check_empty_textures(scn->textures);

    return errs;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection.
scene_intersection intersect_ray(
    const scene* scn, const ray3f& ray, bool find_any) {
    auto iid = 0;
    auto isec = scene_intersection();
    if (!intersect_bvh(
            scn->bvh, ray, find_any, isec.dist, iid, isec.ei, isec.uv))
        return {};
    isec.ist = scn->instances[iid];
    return isec;
}

// Shape element normal.
vec3f eval_elem_norm(const shape* shp, int ei) {
    auto norm = zero3f;
    if (!shp->triangles.empty()) {
        auto t = shp->triangles[ei];
        norm = triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
    } else if (!shp->lines.empty()) {
        auto l = shp->lines[ei];
        norm = line_tangent(shp->pos[l.x], shp->pos[l.y]);
    } else {
        norm = {0, 0, 1};
    }
    return norm;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T eval_elem(const shape* shp, const std::vector<T>& vals, int ei,
    const vec2f& uv, const T& def) {
    if (vals.empty()) return def;
    if (!shp->triangles.empty()) {
        return interpolate_triangle(vals, shp->triangles[ei], uv);
    } else if (!shp->lines.empty()) {
        return interpolate_line(vals, shp->lines[ei], uv.x);
    } else if (!shp->pos.empty()) {
        return vals[ei];
    } else {
        return def;
    }
}

// Shape values interpolated using barycentric coordinates
vec3f eval_pos(const shape* shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->pos, ei, uv, {0, 0, 0});
}
vec3f eval_norm(const shape* shp, int ei, const vec2f& uv) {
    if (shp->norm.empty()) return eval_elem_norm(shp, ei);
    return normalize(eval_elem(shp, shp->norm, ei, uv, {0, 0, 1}));
}
vec2f eval_texcoord(const shape* shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->texcoord, ei, uv, {0, 0});
}
vec4f eval_color(const shape* shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->color, ei, uv, {1, 1, 1, 1});
}
float eval_radius(const shape* shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->radius, ei, uv, 0.0f);
}
vec4f eval_tangsp(const shape* shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->tangsp, ei, uv, {0, 0, 0, 1});
}
vec3f eval_tangsp(
    const shape* shp, int ei, const vec2f& uv, bool& left_handed) {
    auto tangsp = eval_elem(shp, shp->tangsp, ei, uv, {0, 0, 0, 1});
    left_handed = tangsp.w < 0;
    return {tangsp.x, tangsp.y, tangsp.w};
}

// Instance values interpolated using barycentric coordinates.
vec3f eval_pos(const instance* ist, int ei, const vec2f& uv) {
    return transform_point(ist->frame, eval_pos(ist->shp, ei, uv));
}
vec3f eval_norm(const instance* ist, int ei, const vec2f& uv) {
    return transform_direction(ist->frame, eval_norm(ist->shp, ei, uv));
}
vec2f eval_texcoord(const instance* ist, int ei, const vec2f& uv) {
    return eval_texcoord(ist->shp, ei, uv);
}
vec4f eval_color(const instance* ist, int ei, const vec2f& uv) {
    return eval_color(ist->shp, ei, uv);
}
float eval_radius(const instance* ist, int ei, const vec2f& uv) {
    return eval_radius(ist->shp, ei, uv);
}
vec3f eval_tangsp(
    const instance* ist, int ei, const vec2f& uv, bool& left_handed) {
    return transform_direction(
        ist->frame, eval_tangsp(ist->shp, ei, uv, left_handed));
}
// Instance element values.
vec3f eval_elem_norm(const instance* ist, int ei) {
    return transform_direction(ist->frame, eval_elem_norm(ist->shp, ei));
}
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const instance* ist, int ei, const vec2f& uv, const vec3f& o) {
    if (!ist->shp->triangles.empty()) {
        auto n = eval_norm(ist, ei, uv);
        if (ist->mat->norm_txt.txt) {
            auto texcoord = eval_texcoord(ist, ei, uv);
            auto left_handed = false;
            auto tang = eval_tangsp(ist, ei, uv, left_handed);
            auto txt = xyz(eval_texture(ist->mat->norm_txt, texcoord, false));
            txt = txt * 2 - vec3f{1, 1, 1};
            txt.y = -txt.y;
            auto fp = make_frame_fromzx({0, 0, 0}, n, tang);
            if (left_handed) fp.y = -fp.y;
            n = transform_direction(fp, txt);
        }
        if (dot(n, o) < 0) n = -n;
        return n;
    } else if (!ist->shp->lines.empty()) {
        return orthonormalize(o, eval_norm(ist, ei, uv));
    } else {
        return o;
    }
}

// Environment texture coordinates from the direction.
vec2f eval_texcoord(const environment* env, const vec3f& w) {
    auto wl = transform_direction_inverse(env->frame, w);
    auto uv = vec2f{
        atan2(wl.z, wl.x) / (2 * pi), acos(clamp(wl.y, -1.0f, 1.0f)) / pi};
    if (uv.x < 0) uv.x += 1;
    return uv;
}
// Evaluate the environment direction.
vec3f eval_direction(const environment* env, const vec2f& uv) {
    return transform_direction(
        env->frame, {cos(uv.x * 2 * pi) * sin(uv.y * pi), cos(uv.y * pi),
                        sin(uv.x * 2 * pi) * sin(uv.y * pi)});
}
// Evaluate the environment color.
vec3f eval_environment(const environment* env, const vec3f& w) {
    auto ke = env->ke;
    if (env->ke_txt.txt) {
        auto texcoord = eval_texcoord(env, w);
        ke *= xyz(eval_texture(env->ke_txt, texcoord));
    }
    return ke;
}

// Evaluate a texture
vec4f eval_texture(const texture_info& info, const vec2f& texcoord, bool srgb,
    const vec4f& def) {
    auto txt = info.txt;
    if (!txt || (txt->hdr.empty() && txt->ldr.empty())) return def;

    auto lookup = [&def, &txt, &srgb](int i, int j) {
        auto idx = i + j * txt->width;
        if (!txt->ldr.empty())
            return (srgb) ? srgb_to_linear(txt->ldr[idx]) :
                            byte_to_float(txt->ldr[idx]);
        else if (!txt->hdr.empty())
            return txt->hdr[idx];
        else
            return def;
    };

    // get image width/height
    auto w = txt->width, h = txt->height;

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (!info.wrap_s) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * w;
    } else {
        s = std::fmod(texcoord.x, 1.0f) * w;
        if (s < 0) s += w;
    }
    if (!info.wrap_t) {
        t = clamp(texcoord.y, 0.0f, 1.0f) * h;
    } else {
        t = std::fmod(texcoord.y, 1.0f) * h;
        if (t < 0) t += h;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, w - 1), j = clamp((int)t, 0, h - 1);
    auto ii = (i + 1) % w, jj = (j + 1) % h;
    auto u = s - i, v = t - j;

    // nearest lookup
    if (!info.linear) return lookup(i, j);

    // handle interpolation
    return lookup(i, j) * (1 - u) * (1 - v) + lookup(i, jj) * (1 - u) * v +
           lookup(ii, j) * u * (1 - v) + lookup(ii, jj) * u * v;
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv) {
    auto h = 2 * tan(cam->yfov / 2);
    auto w = h * cam->aspect;
    auto o = vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    auto q = vec3f{w * cam->focus * (uv.x - 0.5f),
        h * cam->focus * (uv.y - 0.5f), -cam->focus};
    return make_ray(transform_point(cam->frame, o),
        transform_direction(cam->frame, normalize(q - o)));
}

// Generates a ray from a camera for pixel coordinates `ij`, the
// resolution `res`, the sub-pixel coordinates `puv` and the lens
// coordinates `luv` and the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv) {
    auto uv =
        vec2f{(ij.x + puv.x) / (cam->aspect * res), 1 - (ij.y - puv.y) / res};
    return eval_camera_ray(cam, uv, luv);
}

// Evaluates material parameters.
vec3f eval_emission(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->ke * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->ke_txt, eval_texcoord(ist, ei, uv)));
}
vec3f eval_diffuse(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).x;
        return kb * (1 - km);
    }
}
vec3f eval_specular(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->ks * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).x;
        return kb * km + vec3f{0.04f, 0.04f, 0.04f} * (1 - km);
    }
}
float eval_roughness(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    if (!ist->mat->base_metallic) {
        auto rs = ist->mat->rs;
        auto txt_w =
            eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).w;
        return 1 - ((1 - rs) * txt_w);
    } else {
        return ist->mat->rs *
               eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).w;
    }
    // rs = clamp(rs, 0.03f * 0.03f, 1.0f);
}
vec3f eval_transmission(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->kt * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->kt_txt, eval_texcoord(ist, ei, uv)));
}
float eval_opacity(const instance* ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    return ist->mat->op * eval_color(ist->shp, ei, uv).w *
           eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)).w;
}

// Evaluates the brdf at a location.
brdf eval_brdf(const instance* ist, int ei, const vec2f& uv) {
    auto f = brdf();
    f.kd = eval_diffuse(ist, ei, uv);
    f.ks = eval_specular(ist, ei, uv);
    f.kt = eval_transmission(ist, ei, uv);
    f.rs = eval_roughness(ist, ei, uv);
    if(f.kd != zero3f) {
        f.rs = clamp(f.rs, 0.03f * 0.03f, 1.0f);
    } else if(f.rs <= 0.03f * 0.03f) f.rs = 0;
    return f;
}
bool is_delta_brdf(const brdf& f) {
    return f.rs == 0 && f.kd == zero3f;
}

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const shape* shp, float re, const vec2f& ruv) {
    // TODO: implement sampling without cdf
    if (shp->elem_cdf.empty()) return {};
    if (!shp->triangles.empty()) {
        return sample_triangles(shp->elem_cdf, re, ruv);
    } else if (!shp->lines.empty()) {
        return {sample_lines(shp->elem_cdf, re, ruv.x).first, ruv};
    } else if (!shp->pos.empty()) {
        return {sample_points(shp->elem_cdf, re), ruv};
    } else {
        return {0, zero2f};
    }
}

// Merge scene into one another
void merge_into(scene* merge_into, scene* merge_from) {
    auto merge = [](auto& v1, auto& v2) {
        v1.insert(v1.end(), v2.begin(), v2.end());
        v2.clear();
    };
    merge(merge_into->cameras, merge_from->cameras);
    merge(merge_into->textures, merge_from->textures);
    merge(merge_into->materials, merge_from->materials);
    merge(merge_into->shapes, merge_from->shapes);
    merge(merge_into->environments, merge_from->environments);
    merge(merge_into->nodes, merge_from->nodes);
    merge(merge_into->animations, merge_from->animations);
}

void print_stats(scene* scn) {
    uint64_t num_cameras = 0;
    uint64_t num_shape_groups = 0;
    uint64_t num_shapes = 0;
    uint64_t num_instances = 0;
    uint64_t num_materials = 0;
    uint64_t num_textures = 0;
    uint64_t num_environments = 0;
    uint64_t num_nodes = 0;
    uint64_t num_animations = 0;

    uint64_t elem_lines = 0;
    uint64_t elem_triangles = 0;
    uint64_t vert_pos = 0;
    uint64_t vert_norm = 0;
    uint64_t vert_texcoord = 0;
    uint64_t vert_color = 0;
    uint64_t vert_radius = 0;
    uint64_t vert_tangsp = 0;

    uint64_t texel_ldrs = 0;
    uint64_t texel_hdrs = 0;

    uint64_t memory_ldrs = 0;
    uint64_t memory_hdrs = 0;
    uint64_t memory_elems = 0;
    uint64_t memory_verts = 0;

    update_bbox(scn, true);
    auto bbox = scn->bbox;

    num_cameras = scn->cameras.size();
    num_shapes = scn->shapes.size();
    num_materials = scn->materials.size();
    num_textures = scn->textures.size();
    num_environments = scn->environments.size();
    num_instances = scn->instances.size();
    num_nodes = scn->nodes.size();
    num_animations = scn->animations.size();

    for (auto shp : scn->shapes) {
        elem_lines += shp->lines.size();
        elem_triangles += shp->triangles.size();
        vert_pos += shp->pos.size();
        vert_norm += shp->norm.size();
        vert_texcoord += shp->texcoord.size();
        vert_color += shp->color.size();
        vert_radius += shp->radius.size();
        vert_tangsp += shp->tangsp.size();
    }
    memory_elems = elem_lines * sizeof(vec2i) + elem_triangles * sizeof(vec3i);
    memory_verts = vert_pos * sizeof(vec3f) + vert_norm * sizeof(vec3f) +
                   vert_texcoord * sizeof(vec3f) + vert_color * sizeof(vec4f) +
                   vert_tangsp * sizeof(vec4f) + vert_radius * sizeof(float);

    for (auto txt : scn->textures) {
        texel_ldrs = (txt->ldr.empty()) ? 0 : (txt->width * txt->height);
        texel_hdrs = (txt->hdr.empty()) ? 0 : (txt->width * txt->height);
    }
    memory_ldrs = texel_ldrs * sizeof(vec4b);
    memory_hdrs = texel_hdrs * sizeof(vec4f);

    println("num_cameras: {}", num_cameras);
    println("num_shape_groups: {}", num_shape_groups);
    println("num_shapes: {}", num_shapes);
    println("num_instances: {}", num_instances);
    println("num_materials: {}", num_materials);
    println("num_textures: {}", num_textures);
    println("num_environments: {}", num_environments);
    println("num_nodes: {}", num_nodes);
    println("num_animations: {}", num_animations);
    println("elem_lines: {}", elem_lines);
    println("elem_triangles: {}", elem_triangles);
    println("vert_pos: {}", vert_pos);
    println("vert_norm: {}", vert_norm);
    println("vert_texcoord: {}", vert_texcoord);
    println("vert_color: {}", vert_color);
    println("vert_radius: {}", vert_radius);
    println("vert_tangsp: {}", vert_tangsp);
    println("texel_ldrs: {}", texel_ldrs);
    println("texel_hdrs: {}", texel_hdrs);
    println("memory_ldrs: {}", memory_ldrs);
    println("memory_hdrs: {}", memory_hdrs);
    println("memory_elems: {}", memory_elems);
    println("memory_verts: {}", memory_verts);
    println("bbox_scn: {} {}", bbox.min, bbox.max);
    println("bbox_min   : {}", bbox.min);
    println("bbox_max   : {}", bbox.max);
    println("bbox_size  : {}", bbox.max - bbox.min);
    println("bbox_center: {}", (bbox.max + bbox.min) / 2);
}

#ifdef YGL_OBJ

// Flattens an scene
scene* obj_to_scene(const obj_scene* obj) {
    // clear scene
    auto scn = new scene();

    struct obj_vertex_hash {
        std::hash<int> Th;
        size_t operator()(const obj_vertex& vv) const {
            auto v = (const int*)&vv;
            size_t h = 0;
            for (auto i = 0; i < sizeof(obj_vertex) / sizeof(int); i++) {
                // embads hash_combine below
                h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
            }
            return h;
        }
    };

    // convert textures
    auto tmap = std::unordered_map<std::string, texture*>{{"", nullptr}};
    auto add_texture = [&tmap, scn](const std::string& path) {
        if (tmap.find(path) != tmap.end()) return tmap.at(path);
        auto txt = new texture();
        txt->name = path;
        txt->path = path;
        scn->textures.push_back(txt);
        tmap[path] = txt;
        return txt;
    };

    auto make_texture_info = [&add_texture](const obj_texture_info& oinfo) {
        auto info = texture_info();
        if (oinfo.path == "") return info;
        info.txt = add_texture(oinfo.path);
        info.wrap_s = !oinfo.clamp;
        info.wrap_t = !oinfo.clamp;
        info.scale = oinfo.scale;
        return info;
    };

    // convert materials and build textures
    auto mmap = std::unordered_map<std::string, material*>{{"", nullptr}};
    for (auto omat : obj->materials) {
        auto mat = new material();
        mat->name = omat->name;
        mat->ke = omat->ke;
        mat->kd = omat->kd;
        mat->ks = omat->ks;
        mat->kt = omat->kt;
        if (omat->ns >= 1e6f)
            mat->rs = 0;
        else if (omat->ns < 1)
            mat->rs = 1;
        else
            mat->rs = pow(2 / (omat->ns + 2), 1 / 4.0f);
        mat->op = omat->op;
        mat->ke_txt = make_texture_info(omat->ke_txt);
        mat->kd_txt = make_texture_info(omat->kd_txt);
        mat->ks_txt = make_texture_info(omat->ks_txt);
        mat->kt_txt = make_texture_info(omat->kt_txt);
        mat->norm_txt = make_texture_info(omat->norm_txt);
        mat->bump_txt = make_texture_info(omat->bump_txt);
        mat->disp_txt = make_texture_info(omat->disp_txt);
        scn->materials.push_back(mat);
        mmap[mat->name] = mat;
    }

    // convert meshes
    auto omap = std::unordered_map<std::string,
        std::vector<std::pair<shape*, material*>>>{{"", {}}};
    for (auto omsh : obj->objects) {
        if (omsh->verts.empty()) continue;
        if (omsh->elems.empty()) continue;
        for (auto gid = 0; gid < omsh->groups.size(); gid++) {
            auto ogrp = omsh->groups.at(gid);
            auto shp = new shape();
            shp->name =
                omsh->name + ((gid) ? std::to_string(gid) : std::string());

            // insert all vertices
            std::unordered_map<obj_vertex, int, obj_vertex_hash> vert_map;
            std::vector<int> vert_ids;
            for (auto& vert : omsh->verts) {
                if (vert_map.find(vert) == vert_map.end()) {
                    auto s = (int)vert_map.size();
                    vert_map[vert] = s;
                }
                vert_ids.push_back(vert_map.at(vert));
            }

            // convert elements
            for (auto& elem : omsh->elems) {
                if (elem.groupid != gid) continue;
                switch (elem.type) {
                    case obj_element_type::point: {
                        log_warning("points not supported");
                    } break;
                    case obj_element_type::line: {
                        for (auto i = elem.start;
                             i < elem.start + elem.size - 1; i++) {
                            shp->lines.push_back(
                                {vert_ids[i], vert_ids[i + 1]});
                        }
                    } break;
                    case obj_element_type::face: {
                        for (auto i = elem.start + 2;
                             i < elem.start + elem.size; i++) {
                            shp->triangles.push_back({vert_ids[elem.start],
                                vert_ids[i - 1], vert_ids[i]});
                        }
                    } break;
                    case obj_element_type::bezier: {
                        log_warning("beziers not supported");
                    } break;
                }
            }

            // copy vertex data
            auto v = omsh->verts[0];
            if (v.pos >= 0) shp->pos.resize(vert_map.size());
            if (v.texcoord >= 0) shp->texcoord.resize(vert_map.size());
            if (v.norm >= 0) shp->norm.resize(vert_map.size());
            if (v.color >= 0) shp->color.resize(vert_map.size());
            if (v.radius >= 0) shp->radius.resize(vert_map.size());
            for (auto& kv : vert_map) {
                auto idx = kv.second;
                auto vert = kv.first;
                if (v.pos >= 0 && vert.pos >= 0)
                    shp->pos[idx] = obj->pos[vert.pos];
                if (v.texcoord >= 0 && vert.texcoord >= 0)
                    shp->texcoord[idx] = obj->texcoord[vert.texcoord];
                if (v.norm >= 0 && vert.norm >= 0)
                    shp->norm[idx] = obj->norm[vert.norm];
                if (v.color >= 0 && vert.color >= 0)
                    shp->color[idx] = obj->color[vert.color];
                if (v.radius >= 0 && vert.radius >= 0)
                    shp->radius[idx] = obj->radius[vert.radius];
            }
            scn->shapes.push_back(shp);
            omap[omsh->name].push_back({shp, mmap[ogrp.matname]});
        }
    }

    // convert cameras
    auto cmap = std::unordered_map<std::string, camera*>{{"", nullptr}};
    for (auto ocam : obj->cameras) {
        auto cam = new camera();
        cam->name = ocam->name;
        cam->ortho = ocam->ortho;
        cam->yfov = ocam->yfov;
        cam->aspect = ocam->aspect;
        cam->aperture = ocam->aperture;
        cam->focus = ocam->focus;
        cam->frame = ocam->frame;
        scn->cameras.push_back(cam);
        cmap[cam->name] = cam;
    }

    // convert envs
    std::unordered_set<material*> env_mat;
    auto emap = std::unordered_map<std::string, environment*>{{"", nullptr}};
    for (auto oenv : obj->environments) {
        auto env = new environment();
        env->name = oenv->name;
        env->ke = oenv->ke;
        env->ke_txt = make_texture_info(oenv->ke_txt);
        env->frame = oenv->frame;
        scn->environments.push_back(env);
        emap[env->name] = env;
    }

    // convert nodes
    if (!obj->nodes.empty()) {
        for (auto onde : obj->nodes) {
            auto nde = new node();
            nde->name = onde->name;
            nde->cam = cmap.at(onde->camname);
            nde->env = emap.at(onde->envname);
            nde->translation = onde->translation;
            nde->rotation = onde->rotation;
            nde->scale = onde->scale;
            nde->frame = onde->frame;
            scn->nodes.push_back(nde);
        }

        // set up parent pointers
        for (auto nid = 0; nid < obj->nodes.size(); nid++) {
            auto onde = obj->nodes[nid];
            if (onde->parent.empty()) continue;
            auto nde = scn->nodes[nid];
            for (auto parent : scn->nodes) {
                if (parent->name == onde->parent) {
                    nde->parent = parent;
                    break;
                }
            }
        }

        // set up instances
        for (auto nid = 0; nid < obj->nodes.size(); nid++) {
            auto onde = obj->nodes[nid];
            if (onde->objname.empty()) continue;
            auto nde = scn->nodes[nid];
            auto& shps = omap.at(onde->objname);
            if (shps.empty()) continue;
            if (shps.size() == 1) {
                nde->ist = new instance();
                nde->ist->name = nde->name;
                nde->ist->shp = shps[0].first;
                nde->ist->mat = shps[0].second;
                scn->instances.push_back(nde->ist);
            } else {
                for (auto shp : shps) {
                    auto child = new node();
                    child->name = nde->name + "_" + shp.first->name;
                    child->parent = nde;
                    child->ist = new instance();
                    child->ist->name = child->name;
                    child->ist->shp = shp.first;
                    child->ist->mat = shp.second;
                    scn->instances.push_back(child->ist);
                }
            }
        }
    } else {
        for (auto& shps : omap) {
            for (auto shp : shps.second) {
                auto ist = new instance();
                ist->name = shp.first->name;
                ist->shp = shp.first;
                ist->mat = shp.second;
                scn->instances.push_back(ist);
            }
        }
    }

    // update transforms
    update_transforms(scn);

    // done
    return scn;
}

// Save an scene
obj_scene* scene_to_obj(const scene* scn, bool preserve_instances) {
    auto obj = new obj_scene();

    auto make_texture_info = [](const texture_info& info, bool bump = false) {
        auto oinfo = obj_texture_info();
        if (!info.txt) return oinfo;
        oinfo.path = info.txt->path;
        oinfo.clamp = !info.wrap_s && !info.wrap_t;
        if (bump) oinfo.scale = info.scale;
        return oinfo;
    };

    // convert materials
    for (auto mat : scn->materials) {
        auto omat = new obj_material();
        omat->name = mat->name;
        omat->ke = {mat->ke.x, mat->ke.y, mat->ke.z};
        omat->ke_txt = make_texture_info(mat->ke_txt);
        if (!mat->base_metallic) {
            omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
            omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
            omat->kt = {mat->kt.x, mat->kt.y, mat->kt.z};
            if (mat->rs <= 0)
                omat->ns = 1e6f;
            else if (mat->rs >= 1)
                omat->ns = 0;
            else
                omat->ns = 2 / pow(mat->rs, 4.0f) - 2;
            omat->op = mat->op;
            omat->kd_txt = make_texture_info(mat->kd_txt);
            omat->ks_txt = make_texture_info(mat->ks_txt);
            omat->kt_txt = make_texture_info(mat->kt_txt);
        } else {
            if (mat->rs >= 1 && mat->ks.x == 0) {
                omat->kd = mat->kd;
                omat->ks = {0, 0, 0};
                omat->ns = 1;
            } else {
                auto kd = mat->kd * (1 - 0.04f) * (1 - mat->ks.x);
                auto ks = mat->kd * mat->ks.x +
                          vec3f{0.04f, 0.04f, 0.04f} * (1 - mat->ks.x);
                omat->kd = {kd.x, kd.y, kd.z};
                omat->ks = {ks.x, ks.y, ks.z};
                if (mat->rs <= 0)
                    omat->ns = 1e6f;
                else if (mat->rs >= 1)
                    omat->ns = 0;
                else
                    omat->ns = 2 / pow(mat->rs, 4.0f) - 2;
            }
            omat->op = mat->op;
            if (mat->ks.x < 0.5f) {
                omat->kd_txt = make_texture_info(mat->kd_txt);
            } else {
                omat->ks_txt = make_texture_info(mat->ks_txt);
            }
        }
        omat->bump_txt = make_texture_info(mat->bump_txt, true);
        omat->disp_txt = make_texture_info(mat->disp_txt, true);
        omat->norm_txt = make_texture_info(mat->norm_txt, true);
        if (mat->op < 1 || mat->kt != zero3f) {
            omat->illum = 4;
        } else {
            omat->illum = 2;
        }
        obj->materials.push_back(omat);
    }

    // add elem
    auto add_elem = [](auto& shp, auto& oobj, auto etype, auto esize, auto ei) {
        auto elem = obj_element();
        elem.start = (uint32_t)oobj->verts.size();
        elem.type = etype;
        elem.size = (uint16_t)esize;
        elem.groupid = 0;
        oobj->elems.push_back(elem);
    };
    // add vertex
    auto add_vert = [](auto& shp, auto& oobj, auto& offset, int vid) {
        auto vert = obj_vertex{-1, -1, -1, -1, -1};
        if (!shp->pos.empty()) vert.pos = offset.pos + vid;
        if (!shp->texcoord.empty()) vert.texcoord = offset.texcoord + vid;
        if (!shp->norm.empty()) vert.norm = offset.norm + vid;
        if (!shp->color.empty()) vert.color = offset.color + vid;
        if (!shp->radius.empty()) vert.radius = offset.radius + vid;
        oobj->verts.push_back(vert);
    };

    // flatten instances if necessary
    auto flatten_instances =
        !preserve_instances && scn->nodes.empty() && !scn->instances.empty();
    auto shapes = std::vector<std::pair<shape*, frame3f>>();
    if (flatten_instances) {
        for (auto ist : scn->instances)
            shapes.push_back({ist->shp, ist->frame});
    } else {
        for (auto shp : scn->shapes) shapes.push_back({shp, identity_frame3f});
    }
    auto shape_mats = std::map<shape*, material*>();
    for (auto shp : scn->shapes) shape_mats[shp] = nullptr;
    for (auto ist : scn->instances) {
        if (shape_mats.at(ist->shp) && shape_mats.at(ist->shp) != ist->mat)
            log_error("shapes can only have one material associated");
        else
            shape_mats[ist->shp] = ist->mat;
    }

    // convert shapes
    for (auto& shp_frame : shapes) {
        auto shp = shp_frame.first;
        auto frame = shp_frame.second;
        auto mat = shape_mats.at(shp);
        auto oobj = new obj_object();
        oobj->name = shp->name;
        auto offset = obj_vertex{(int)obj->pos.size(),
            (int)obj->texcoord.size(), (int)obj->norm.size(),
            (int)obj->color.size(), (int)obj->radius.size()};
        if (frame != identity_frame3f) {
            for (auto& v : shp->pos)
                obj->pos.push_back(transform_point(frame, v));
            for (auto& v : shp->norm)
                obj->norm.push_back(transform_direction(frame, v));
        } else {
            for (auto& v : shp->pos) obj->pos.push_back(v);
            for (auto& v : shp->norm) obj->norm.push_back(v);
        }
        for (auto& v : shp->texcoord) obj->texcoord.push_back(v);
        for (auto& v : shp->color) obj->color.push_back({v.x, v.y, v.z, v.w});
        for (auto& v : shp->radius) obj->radius.push_back(v);
        oobj->groups.push_back(
            {"", (mat) ? mat->name : std::string(), shp->norm.empty()});
        for (auto ei = 0; ei < shp->lines.size(); ei++) {
            auto l = shp->lines[ei];
            add_elem(shp, oobj, obj_element_type::line, 2, ei);
            for (auto vid : {l.x, l.y}) add_vert(shp, oobj, offset, vid);
        }
        for (auto ei = 0; ei < shp->triangles.size(); ei++) {
            auto t = shp->triangles[ei];
            add_elem(shp, oobj, obj_element_type::face, 3, ei);
            for (auto vid : {t.x, t.y, t.z}) add_vert(shp, oobj, offset, vid);
        }
        obj->objects.push_back(oobj);
    }

    // convert cameras
    for (auto cam : scn->cameras) {
        auto ocam = new obj_camera();
        ocam->name = cam->name;
        ocam->ortho = cam->ortho;
        ocam->yfov = cam->yfov;
        ocam->aspect = cam->aspect;
        ocam->focus = cam->focus;
        ocam->aperture = cam->aperture;
        ocam->frame = cam->frame;
        obj->cameras.push_back(ocam);
    }

    // convert envs
    for (auto env : scn->environments) {
        auto oenv = new obj_environment();
        oenv->name = env->name;
        oenv->ke = env->ke;
        oenv->ke_txt = make_texture_info(env->ke_txt);
        oenv->frame = env->frame;
        obj->environments.push_back(oenv);
    }

    // convert hierarchy
    if (!scn->nodes.empty()) {
        for (auto nde : scn->nodes) {
            auto onde = new obj_node();
            onde->name = nde->name;
            if (nde->cam) onde->camname = nde->cam->name;
            if (nde->ist) onde->objname = nde->ist->shp->name;
            if (nde->env) onde->envname = nde->env->name;
            onde->frame = nde->frame;
            onde->translation = nde->translation;
            onde->rotation = nde->rotation;
            onde->scale = nde->scale;
            obj->nodes.push_back(onde);
        }

        // parent
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto onde = obj->nodes.at(idx);
            onde->parent = nde->parent->name;
        }
    } else if (!flatten_instances) {
        for (auto ist : scn->instances) {
            auto onde = new obj_node();
            onde->name = ist->name;
            onde->objname = ist->shp->name;
            onde->frame = ist->frame;
            obj->nodes.push_back(onde);
        }
    }

    return obj;
}

#endif

#ifdef YGL_GLTF

// Flattens a gltf file into a flattened asset.
scene* gltf_to_scene(const glTF* gltf) {
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };

    // clear asset
    auto scn = new scene();

    // convert images
    for (auto gtxt : gltf->images) {
        auto txt = new texture();
        txt->name = gtxt->name;
        txt->path = (startswith(gtxt->uri, "data:")) ?
                        std::string("[glTF inline]") :
                        gtxt->uri;
        scn->textures.push_back(txt);
    }

    // add a texture
    auto make_texture_info = [gltf, scn](const glTFTextureInfo* ginfo,
                                 bool normal = false, bool occlusion = false) {
        auto info = texture_info();
        if (!ginfo) return info;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return info;
        info.txt = scn->textures.at((int)gtxt->source);
        if (!info.txt) return info;
        auto gsmp = gltf->get(gtxt->sampler);
        if (gsmp) {
            info.linear = gsmp->magFilter != glTFSamplerMagFilter::Nearest;
            info.mipmap = gsmp->minFilter != glTFSamplerMinFilter::Linear &&
                          gsmp->minFilter != glTFSamplerMinFilter::Nearest;
            info.wrap_s = gsmp->wrapS != glTFSamplerWrapS::ClampToEdge;
            info.wrap_t = gsmp->wrapT != glTFSamplerWrapT::ClampToEdge;
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            info.scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            info.scale = ninfo->strength;
        }
        return info;
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new material();
        mat->name = gmat->name;
        mat->ke = gmat->emissiveFactor;
        mat->ke_txt = make_texture_info(gmat->emissiveTexture);
        if (gmat->pbrSpecularGlossiness) {
            mat->base_metallic = false;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor.x, gsg->diffuseFactor.y,
                gsg->diffuseFactor.z};
            mat->op = gsg->diffuseFactor.w;
            mat->ks = gsg->specularFactor;
            mat->rs = 1 - gsg->glossinessFactor;
            mat->kd_txt = make_texture_info(gsg->diffuseTexture);
            mat->ks_txt = make_texture_info(gsg->specularGlossinessTexture);
        } else if (gmat->pbrMetallicRoughness) {
            mat->base_metallic = true;
            auto gmr = gmat->pbrMetallicRoughness;
            mat->kd = {gmr->baseColorFactor.x, gmr->baseColorFactor.y,
                gmr->baseColorFactor.z};
            mat->op = gmr->baseColorFactor.w;
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = make_texture_info(gmr->baseColorTexture);
            mat->ks_txt = make_texture_info(gmr->metallicRoughnessTexture);
        }
        mat->norm_txt = make_texture_info(gmat->normalTexture, true, false);
        scn->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<std::pair<shape*, material*>>>();
    for (auto gmesh : gltf->meshes) {
        meshes.push_back({});
        auto sid = 0;
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
            shp->name =
                gmesh->name + ((sid) ? std::to_string(sid) : std::string());
            sid++;
            for (auto gattr : gprim->attributes) {
                auto semantic = gattr.first;
                auto vals = accessor_view(gltf, gltf->get(gattr.second));
                if (semantic == "POSITION") {
                    shp->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->pos.push_back(vals.getv3f(i));
                } else if (semantic == "NORMAL") {
                    shp->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->norm.push_back(vals.getv3f(i));
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    shp->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->texcoord.push_back(vals.getv2f(i));
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    shp->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->color.push_back(vals.getv4f(i, {0, 0, 0, 1}));
                } else if (semantic == "TANGENT") {
                    shp->tangsp.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->tangsp.push_back(vals.getv4f(i));
                } else if (semantic == "RADIUS") {
                    shp->radius.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->radius.push_back(vals.get(i, 0));
                } else {
                    // ignore
                }
            }
            // indices
            if (!gprim->indices) {
                if (gprim->mode == glTFMeshPrimitiveMode::Triangles) {
                    shp->triangles.reserve(shp->pos.size() / 3);
                    for (auto i = 0; i < shp->pos.size() / 3; i++)
                        shp->triangles.push_back(
                            {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                } else if (gprim->mode == glTFMeshPrimitiveMode::TriangleFan) {
                    shp->triangles.reserve(shp->pos.size() - 2);
                    for (auto i = 2; i < shp->pos.size(); i++)
                        shp->triangles.push_back({0, i - 1, i});
                } else if (gprim->mode ==
                           glTFMeshPrimitiveMode::TriangleStrip) {
                    shp->triangles.reserve(shp->pos.size() - 2);
                    for (auto i = 2; i < shp->pos.size(); i++)
                        shp->triangles.push_back({i - 2, i - 1, i});
                } else if (gprim->mode == glTFMeshPrimitiveMode::Lines) {
                    shp->lines.reserve(shp->pos.size() / 2);
                    for (auto i = 0; i < shp->pos.size() / 2; i++)
                        shp->lines.push_back({i * 2 + 0, i * 2 + 1});
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineLoop) {
                    shp->lines.reserve(shp->pos.size());
                    for (auto i = 1; i < shp->pos.size(); i++)
                        shp->lines.push_back({i - 1, i});
                    shp->lines.back() = {(int)shp->pos.size() - 1, 0};
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineStrip) {
                    shp->lines.reserve(shp->pos.size() - 1);
                    for (auto i = 1; i < shp->pos.size(); i++)
                        shp->lines.push_back({i - 1, i});
                } else if (gprim->mode == glTFMeshPrimitiveMode::NotSet ||
                           gprim->mode == glTFMeshPrimitiveMode::Points) {
                    log_warning("points not supported");
                } else {
                    throw std::runtime_error("unknown primitive type");
                }
            } else {
                auto indices = accessor_view(gltf, gltf->get(gprim->indices));
                if (gprim->mode == glTFMeshPrimitiveMode::Triangles) {
                    shp->triangles.reserve(indices.size());
                    for (auto i = 0; i < indices.size() / 3; i++)
                        shp->triangles.push_back({indices.geti(i * 3 + 0),
                            indices.geti(i * 3 + 1), indices.geti(i * 3 + 2)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::TriangleFan) {
                    shp->triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shp->triangles.push_back({indices.geti(0),
                            indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode ==
                           glTFMeshPrimitiveMode::TriangleStrip) {
                    shp->triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shp->triangles.push_back({indices.geti(i - 2),
                            indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::Lines) {
                    shp->lines.reserve(indices.size() / 2);
                    for (auto i = 0; i < indices.size() / 2; i++)
                        shp->lines.push_back(
                            {indices.geti(i * 2 + 0), indices.geti(i * 2 + 1)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineLoop) {
                    shp->lines.reserve(indices.size());
                    for (auto i = 1; i < indices.size(); i++)
                        shp->lines.push_back(
                            {indices.geti(i - 1), indices.geti(i)});
                    shp->lines.back() = {
                        indices.geti(indices.size() - 1), indices.geti(0)};
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineStrip) {
                    shp->lines.reserve(indices.size() - 1);
                    for (auto i = 1; i < indices.size(); i++)
                        shp->lines.push_back(
                            {indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::NotSet ||
                           gprim->mode == glTFMeshPrimitiveMode::Points) {
                    log_warning("points not supported");
                } else {
                    throw std::runtime_error("unknown primitive type");
                }
            }
            meshes.back().push_back(
                {shp, scn->materials[(int)gprim->material]});
            scn->shapes.push_back(shp);
        }
    }

    // convert cameras
    for (auto gcam : gltf->cameras) {
        auto cam = new camera();
        cam->name = gcam->name;
        cam->ortho = gcam->type == glTFCameraType::Orthographic;
        if (cam->ortho) {
            auto ortho = gcam->orthographic;
            cam->yfov = ortho->ymag;
            cam->aspect = ortho->xmag / ortho->ymag;
            cam->near = ortho->znear;
            cam->far = ortho->zfar;
        } else {
            auto persp = gcam->perspective;
            cam->yfov = persp->yfov;
            cam->aspect = persp->aspectRatio;
            if (!cam->aspect) cam->aspect = 16.0f / 9.0f;
            cam->near = persp->znear;
            cam->far = persp->zfar;
        }
        scn->cameras.push_back(cam);
    }

    // convert nodes
    for (auto gnde : gltf->nodes) {
        auto nde = new node();
        nde->name = gnde->name;
        if (gnde->camera) nde->cam = scn->cameras[(int)gnde->camera];
        nde->translation = gnde->translation;
        nde->rotation = gnde->rotation;
        nde->scale = gnde->scale;
        nde->frame = mat_to_frame(gnde->matrix);
        scn->nodes.push_back(nde);
    }

    // set up parent pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        auto nde = scn->nodes[nid];
        for (auto cid : gnde->children) scn->nodes[(int)cid]->parent = nde;
    }

    // set up instances
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        if (!gnde->mesh) continue;
        auto nde = scn->nodes[nid];
        auto& shps = meshes.at((int)gnde->mesh);
        if (shps.empty()) continue;
        if (shps.size() == 1) {
            nde->ist = new instance();
            nde->ist->name = nde->name;
            nde->ist->shp = shps[0].first;
            nde->ist->mat = shps[0].second;
            scn->instances.push_back(nde->ist);
        } else {
            for (auto shp : shps) {
                auto child = new node();
                child->name = nde->name + "_" + shp.first->name;
                child->parent = nde;
                child->ist = new instance();
                child->ist->name = child->name;
                child->ist->shp = shp.first;
                child->ist->mat = shp.second;
                scn->instances.push_back(child->ist);
            }
        }
    }

    // keyframe type conversion
    static auto keyframe_types =
        std::unordered_map<glTFAnimationSamplerInterpolation, animation_type>{
            {glTFAnimationSamplerInterpolation::NotSet, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Linear, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Step, animation_type::step},
            {glTFAnimationSamplerInterpolation::CubicSpline,
                animation_type::bezier},
        };

    // convert animations
    for (auto ganm : gltf->animations) {
        auto aid = 0;
        auto sampler_map = std::unordered_map<vec2i, int>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto anm = new animation();
                anm->name = ((ganm->name != "") ? ganm->name : "anim") +
                            std::to_string(aid++);
                anm->group = ganm->name;
                auto input_view =
                    accessor_view(gltf, gltf->get(gsampler->input));
                anm->times.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    anm->times[i] = input_view.get(i);
                anm->type = keyframe_types.at(gsampler->interpolation);
                auto output_view =
                    accessor_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        anm->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->translation.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        anm->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->rotation.push_back(output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        anm->scale.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->scale.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Weights: {
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp = max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = std::vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            anm->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < anm->weights.size(); i++) {
                                anm->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    anm->weights[i][j] = values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        throw std::runtime_error("should not have gotten here");
                    }
                }
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = (int)scn->animations.size();
                scn->animations.push_back(anm);
            }
            scn->animations[sampler_map.at({(int)gchannel->sampler,
                                (int)gchannel->target->path})]
                ->targets.push_back(scn->nodes[(int)gchannel->target->node]);
        }
    }

    // compute transforms
    update_transforms(scn, 0);

    return scn;
}

// Unflattnes gltf
glTF* scene_to_gltf(
    const scene* scn, const std::string& buffer_uri, bool separate_buffers) {
    auto gltf = new glTF();

    // add asset info
    gltf->asset = new glTFAsset();
    gltf->asset->generator = "Yocto/gltf";
    gltf->asset->version = "2.0";

    // convert cameras
    for (auto cam : scn->cameras) {
        auto gcam = new glTFCamera();
        gcam->name = cam->name;
        gcam->type = (cam->ortho) ? glTFCameraType::Orthographic :
                                    glTFCameraType::Perspective;
        if (cam->ortho) {
            auto ortho = new glTFCameraOrthographic();
            ortho->ymag = cam->yfov;
            ortho->xmag = cam->aspect * cam->yfov;
            ortho->znear = cam->near;
            ortho->znear = cam->far;
            gcam->orthographic = ortho;
        } else {
            auto persp = new glTFCameraPerspective();
            persp->yfov = cam->yfov;
            persp->aspectRatio = cam->aspect;
            persp->znear = cam->near;
            persp->zfar = cam->far;
            gcam->perspective = persp;
        }
        gltf->cameras.push_back(gcam);
    }

    // convert images
    for (auto txt : scn->textures) {
        auto gimg = new glTFImage();
        gimg->uri = txt->path;
        gltf->images.push_back(gimg);
    }

    // index of an object
    auto index = [](const auto& vec, auto& val) -> int {
        auto pos = find(vec.begin(), vec.end(), val);
        if (pos == vec.end()) return -1;
        return (int)(pos - vec.begin());
    };

    // add a texture and sampler
    auto add_texture_info = [&gltf, &index, scn](const texture_info& info,
                                bool norm = false, bool occ = false) {
        if (!info.txt) return (glTFTextureInfo*)nullptr;
        auto gtxt = new glTFTexture();
        gtxt->name = info.txt->name;
        gtxt->source = glTFid<glTFImage>(index(scn->textures, info.txt));

        // check if it is default
        auto is_default =
            info.wrap_s && info.wrap_t && info.linear && info.mipmap;
        if (!is_default) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = (info.wrap_s) ? glTFSamplerWrapS::Repeat :
                                          glTFSamplerWrapS::ClampToEdge;
            gsmp->wrapT = (info.wrap_t) ? glTFSamplerWrapT::Repeat :
                                          glTFSamplerWrapT::ClampToEdge;
            gsmp->minFilter = (info.mipmap) ?
                                  glTFSamplerMinFilter::LinearMipmapLinear :
                                  glTFSamplerMinFilter::Nearest;
            gsmp->magFilter = (info.linear) ? glTFSamplerMagFilter::Linear :
                                              glTFSamplerMagFilter::Nearest;
            gtxt->sampler = glTFid<glTFSampler>((int)gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        if (norm) {
            auto ginfo = new glTFMaterialNormalTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->scale = info.scale;
            return (glTFTextureInfo*)ginfo;
        } else if (occ) {
            auto ginfo = new glTFMaterialOcclusionTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->strength = info.scale;
            return (glTFTextureInfo*)ginfo;
        } else {
            auto ginfo = new glTFTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            return ginfo;
        }
    };

    // convert materials
    for (auto mat : scn->materials) {
        auto gmat = new glTFMaterial();
        gmat->name = mat->name;
        gmat->emissiveFactor = mat->ke;
        gmat->emissiveTexture = add_texture_info(mat->ke_txt);
        if (!mat->base_metallic) {
            gmat->pbrSpecularGlossiness =
                new glTFMaterialPbrSpecularGlossiness();
            auto gsg = gmat->pbrSpecularGlossiness;
            gsg->diffuseFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
            gsg->specularFactor = mat->ks;
            gsg->glossinessFactor = 1 - mat->rs;
            gsg->diffuseTexture = add_texture_info(mat->kd_txt);
            gsg->specularGlossinessTexture = add_texture_info(mat->ks_txt);
        } else {
            gmat->pbrMetallicRoughness = new glTFMaterialPbrMetallicRoughness();
            auto gmr = gmat->pbrMetallicRoughness;
            gmr->baseColorFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
            gmr->metallicFactor = mat->ks.x;
            gmr->roughnessFactor = mat->rs;
            gmr->baseColorTexture = add_texture_info(mat->kd_txt);
            gmr->metallicRoughnessTexture = add_texture_info(mat->ks_txt);
        }
        gmat->normalTexture = (glTFMaterialNormalTextureInfo*)(add_texture_info(
            mat->norm_txt, true, false));
        gmat->doubleSided = true;
        gltf->materials.push_back(gmat);
    }

    // add buffer
    auto add_buffer = [&gltf](const std::string& buffer_uri) {
        auto gbuffer = new glTFBuffer();
        gltf->buffers.push_back(gbuffer);
        gbuffer->uri = buffer_uri;
        return gbuffer;
    };

    // init buffers
    auto gbuffer_global = add_buffer(buffer_uri);

    // add an optional buffer
    auto add_opt_buffer = [&gbuffer_global, buffer_uri, &add_buffer,
                              separate_buffers](const std::string& uri) {
        if (separate_buffers && uri != "") {
            return add_buffer(uri);
        } else {
            if (!gbuffer_global) gbuffer_global = add_buffer(buffer_uri);
            return gbuffer_global;
        }
    };

    // attribute handling
    auto add_accessor = [&gltf, &index](glTFBuffer* gbuffer,
                            const std::string& name, glTFAccessorType type,
                            glTFAccessorComponentType ctype, int count,
                            int csize, const void* data, bool save_min_max) {
        gltf->bufferViews.push_back(new glTFBufferView());
        auto bufferView = gltf->bufferViews.back();
        bufferView->buffer = glTFid<glTFBuffer>(index(gltf->buffers, gbuffer));
        bufferView->byteOffset = (int)gbuffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = count * csize;
        gbuffer->data.resize(gbuffer->data.size() + bufferView->byteLength);
        gbuffer->byteLength += bufferView->byteLength;
        auto ptr = gbuffer->data.data() + gbuffer->data.size() -
                   bufferView->byteLength;
        bufferView->target = glTFBufferViewTarget::ArrayBuffer;
        memcpy(ptr, data, bufferView->byteLength);
        gltf->accessors.push_back(new glTFAccessor());
        auto accessor = gltf->accessors.back();
        accessor->bufferView =
            glTFid<glTFBufferView>((int)gltf->bufferViews.size() - 1);
        accessor->byteOffset = 0;
        accessor->componentType = ctype;
        accessor->count = count;
        accessor->type = type;
        if (save_min_max && count &&
            ctype == glTFAccessorComponentType::Float) {
            float dmin[4] = {flt_max, flt_max, flt_max, flt_max};
            float dmax[4] = {flt_min, flt_min, flt_min, flt_min};
            auto d = (float*)data;
            auto nc = 0;
            switch (type) {
                case glTFAccessorType::Scalar: nc = 1; break;
                case glTFAccessorType::Vec2: nc = 2; break;
                case glTFAccessorType::Vec3: nc = 3; break;
                case glTFAccessorType::Vec4: nc = 4; break;
                default: break;
            }
            for (auto i = 0; i < count; i++) {
                for (auto c = 0; c < nc; c++) {
                    dmin[c] = min(dmin[c], d[i * nc + c]);
                    dmax[c] = max(dmax[c], d[i * nc + c]);
                }
            }
            for (auto c = 0; c < nc; c++) {
                accessor->min.push_back(dmin[c]);
                accessor->max.push_back(dmax[c]);
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // instances
    auto shape_mats = std::map<shape*, material*>();
    for (auto shp : scn->shapes) shape_mats[shp] = nullptr;
    for (auto ist : scn->instances) {
        if (shape_mats.at(ist->shp) && shape_mats.at(ist->shp) != ist->mat)
            log_error("shapes can only have one material associated");
        else
            shape_mats[ist->shp] = ist->mat;
    }

    // convert meshes
    for (auto shp : scn->shapes) {
        auto gmesh = new glTFMesh();
        gmesh->name = shp->name;
        auto gbuffer = add_opt_buffer(shp->path);
        auto gprim = new glTFMeshPrimitive();
        gprim->material =
            glTFid<glTFMaterial>(index(scn->materials, shape_mats.at(shp)));
        if (!shp->pos.empty())
            gprim->attributes["POSITION"] =
                add_accessor(gbuffer, shp->name + "_pos",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)shp->pos.size(), sizeof(vec3f), shp->pos.data(), true);
        if (!shp->norm.empty())
            gprim->attributes["NORMAL"] = add_accessor(gbuffer,
                shp->name + "_norm", glTFAccessorType::Vec3,
                glTFAccessorComponentType::Float, (int)shp->norm.size(),
                sizeof(vec3f), shp->norm.data(), false);
        if (!shp->texcoord.empty())
            gprim->attributes["TEXCOORD_0"] = add_accessor(gbuffer,
                shp->name + "_texcoord", glTFAccessorType::Vec2,
                glTFAccessorComponentType::Float, (int)shp->texcoord.size(),
                sizeof(vec2f), shp->texcoord.data(), false);
        if (!shp->color.empty())
            gprim->attributes["COLOR_0"] = add_accessor(gbuffer,
                shp->name + "_color", glTFAccessorType::Vec4,
                glTFAccessorComponentType::Float, (int)shp->color.size(),
                sizeof(vec4f), shp->color.data(), false);
        if (!shp->radius.empty())
            gprim->attributes["RADIUS"] = add_accessor(gbuffer,
                shp->name + "_radius", glTFAccessorType::Scalar,
                glTFAccessorComponentType::Float, (int)shp->radius.size(),
                sizeof(float), shp->radius.data(), false);
        if (!shp->lines.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_lines",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->lines.size() * 2, sizeof(int),
                (int*)shp->lines.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Lines;
        } else if (!shp->triangles.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_triangles",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->triangles.size() * 3, sizeof(int),
                (int*)shp->triangles.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Triangles;
        }
        gmesh->primitives.push_back(gprim);
        gltf->meshes.push_back(gmesh);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // shapes
        for (auto ist : scn->instances) {
            auto gnode = new glTFNode();
            gnode->name = ist->name;
            gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, ist->shp));
            gnode->matrix = frame_to_mat(ist->frame);
            gltf->nodes.push_back(gnode);
        }

        // cameras
        for (auto cam : scn->cameras) {
            auto gnode = new glTFNode();
            gnode->name = cam->name;
            gnode->camera = glTFid<glTFCamera>(index(scn->cameras, cam));
            gnode->matrix = frame_to_mat(cam->frame);
            gltf->nodes.push_back(gnode);
        }

        // scenes
        if (!gltf->nodes.empty()) {
            auto gscene = new glTFScene();
            gscene->name = "scene";
            for (auto i = 0; i < gltf->nodes.size(); i++) {
                gscene->nodes.push_back(glTFid<glTFNode>(i));
            }
            gltf->scenes.push_back(gscene);
            gltf->scene = glTFid<glTFScene>(0);
        }
    } else {
        for (auto nde : scn->nodes) {
            auto gnode = new glTFNode();
            gnode->name = nde->name;
            if (nde->cam) {
                gnode->camera =
                    glTFid<glTFCamera>(index(scn->cameras, nde->cam));
            }
            if (nde->ist) {
                gnode->mesh =
                    glTFid<glTFMesh>(index(scn->shapes, nde->ist->shp));
            }
            gnode->matrix = frame_to_mat(nde->frame);
            gnode->translation = nde->translation;
            gnode->rotation = nde->rotation;
            gnode->scale = nde->scale;
            gltf->nodes.push_back(gnode);
        }

        // children
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto gnde = gltf->nodes.at(index(scn->nodes, nde->parent));
            gnde->children.push_back(glTFid<glTFNode>(idx));
        }

        // root nodes
        auto is_root = std::vector<bool>(gltf->nodes.size(), true);
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            auto gnde = gltf->nodes.at(idx);
            for (auto idx1 = 0; idx1 < gnde->children.size(); idx1++) {
                is_root[(int)gnde->children.at(idx1)] = false;
            }
        }

        // scene with root nodes
        auto gscene = new glTFScene();
        gscene->name = "scene";
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            if (is_root[idx]) gscene->nodes.push_back(glTFid<glTFNode>(idx));
        }
        gltf->scenes.push_back(gscene);
        gltf->scene = glTFid<glTFScene>(0);
    }

    // interpolation map
    static const auto interpolation_map =
        std::map<animation_type, glTFAnimationSamplerInterpolation>{
            {animation_type::step, glTFAnimationSamplerInterpolation::Step},
            {animation_type::linear, glTFAnimationSamplerInterpolation::Linear},
            {animation_type::bezier,
                glTFAnimationSamplerInterpolation::CubicSpline},
        };

    // gruop animations
    struct anim_group {
        std::string path;
        std::string name;
        std::vector<animation*> animations;
    };

    std::map<std::string, anim_group> anim_groups;
    for (auto anm : scn->animations) {
        auto agr = &anim_groups[anm->group];
        if (agr->path == "") agr->path = anm->path;
        agr->animations.push_back(anm);
    }

    // animation
    for (auto& agr_kv : anim_groups) {
        auto agr = &agr_kv.second;
        auto ganm = new glTFAnimation();
        ganm->name = agr->name;
        auto gbuffer = add_opt_buffer(agr->path);
        auto count = 0;
        for (auto anm : scn->animations) {
            auto aid = ganm->name + "_" + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input =
                add_accessor(gbuffer, aid + "_time", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)anm->times.size(),
                    sizeof(float), anm->times.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!anm->translation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->translation.size(), sizeof(vec3f),
                    anm->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!anm->rotation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)anm->rotation.size(), sizeof(vec4f),
                    anm->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!anm->scale.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->scale.size(), sizeof(vec3f), anm->scale.data(),
                    false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!anm->weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(anm->weights.size() * anm->weights[0].size());
                for (auto i = 0; i < anm->weights.size(); i++) {
                    values.insert(values.end(), anm->weights[i].begin(),
                        anm->weights[i].end());
                }
                gsmp->output = add_accessor(gbuffer, aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
                throw std::runtime_error("should not have gotten here");
            }
            gsmp->interpolation = interpolation_map.at(anm->type);
            for (auto target : anm->targets) {
                auto gchan = new glTFAnimationChannel();
                gchan->sampler =
                    glTFid<glTFAnimationSampler>{(int)ganm->samplers.size()};
                gchan->target = new glTFAnimationChannelTarget();
                gchan->target->node =
                    glTFid<glTFNode>{index(scn->nodes, target)};
                gchan->target->path = path;
                ganm->channels.push_back(gchan);
            }
            ganm->samplers.push_back(gsmp);
        }

        gltf->animations.push_back(ganm);
    }

    // done
    return gltf;
}

#endif

// Loads/saves textures
void load_textures(const std::string& filename, scene* scn, bool skip_missing) {
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        // TODO: handle glTF buffer textures
        if (txt->path == "[glTF inline]") {
            log_warning("cannot handle glTF inline images");
            txt->width = 1;
            txt->height = 1;
            txt->ldr = {vec4b{255, 255, 255, 255}};
            continue;
        }
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
#if YGL_IMAGEIO
        if (is_hdr_filename(filename)) {
            txt->hdr = load_image4f(filename, txt->width, txt->height);
        } else {
            txt->ldr = load_image4b(filename, txt->width, txt->height);
        }
#endif
        if (txt->hdr.empty() && txt->ldr.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}
void save_textures(
    const std::string& filename, const scene* scn, bool skip_missing) {
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        if (txt->ldr.empty() && txt->hdr.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
#if YGL_IMAGEIO
        if (!txt->ldr.empty()) {
            ok = save_image4b(filename, txt->width, txt->height, txt->ldr);
        }
        if (!txt->hdr.empty()) {
            ok = save_image4f(filename, txt->width, txt->height, txt->hdr);
        }
#endif
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
}

// Load a scene
scene* load_scene(const std::string& filename, bool load_txts,
    bool split_obj_shapes, bool skip_missing) {
    auto ext = path_extension(filename);
    auto scn = (scene*)nullptr;
    if (ext == ".obj" || ext == ".OBJ") {
#if YGL_OBJ
        auto oscn = load_obj(filename, split_obj_shapes, true, true);
        scn = obj_to_scene(oscn);
        scn->name = path_filename(filename);
        delete oscn;
#else
        throw std::runtime_error("Obj not supported");
#endif
    } else if (ext == ".gltf" || ext == ".GLTF") {
#if YGL_GLTF
        auto gscn = load_gltf(filename, true);
        scn = gltf_to_scene(gscn);
        delete gscn;
#else
        throw std::runtime_error("glTF not supported");
#endif
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
    if (scn->name == "") scn->name = path_filename(filename);
    if (!scn) throw std::runtime_error("could not convert gltf scene");
    if (load_txts) load_textures(filename, scn, skip_missing);
    return scn;
}

// Save a scene
void save_scene(const std::string& filename, const scene* scn, bool save_txt,
    bool preserve_obj_instances, bool gltf_separate_buffers,
    bool skip_missing) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ") {
#if YGL_OBJ
        auto oscn = scene_to_obj(scn, preserve_obj_instances);
        save_obj(filename, oscn, true, false);
        delete oscn;
#else
        throw std::runtime_error("unsupported Obj");
#endif
    } else if (ext == ".gltf" || ext == ".GLTF") {
#if YGL_GLTF
        auto buffer_uri = path_basename(filename) + ".bin";
        auto gscn = scene_to_gltf(scn, buffer_uri, gltf_separate_buffers);
        save_gltf(filename, gscn);
        delete gscn;
#else
        throw std::runtime_error("unsupported glTF");
#endif
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
    if (save_txt) save_textures(filename, scn, skip_missing);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

camera* make_camera(const std::string& name, const vec3f& from, const vec3f& to,
    float yfov, float aspect) {
    auto cam = new camera();
    cam->name = name;
    cam->frame = lookat_frame(from, to, vec3f{0, 1, 0});
    cam->yfov = yfov;
    cam->aspect = aspect;
    cam->near = 0.01f;
    cam->far = 10000;
    cam->aperture = 0;
    cam->focus = length(from - to);
    return cam;
};

material* make_material(
    const std::string& name, const vec3f& kd, const vec3f& ks, float rs) {
    auto mat = new material();
    mat->name = name;
    mat->kd = kd;
    mat->ks = ks;
    mat->rs = rs;
    return mat;
}

texture* make_texture(const std::string& name, const std::string& path,
    int width, int height, const std::vector<vec4b>& ldr,
    const std::vector<vec4f>& hdr) {
    auto txt = new texture();
    txt->name = name;
    txt->path = path;
    txt->width = width;
    txt->height = height;
    txt->ldr = ldr;
    txt->hdr = hdr;
    return txt;
}

instance* make_instance(
    const std::string& name, shape* shp, material* mat, const frame3f& frame) {
    auto ist = new instance();
    ist->name = name;
    ist->shp = shp;
    ist->mat = mat;
    ist->frame = frame;
    return ist;
}

node* make_node(const std::string& name, camera* cam, instance* ist,
    environment* env, const frame3f& frame) {
    auto nde = new node();
    nde->name = name;
    nde->cam = cam;
    nde->ist = ist;
    nde->env = env;
    nde->frame = frame;
    return nde;
}

environment* make_environment(const std::string& name, const vec3f& ke,
    texture* ke_txt, const frame3f& frame) {
    auto env = new environment();
    env->name = name;
    env->ke = ke;
    env->ke_txt.txt = ke_txt;
    env->frame = frame;
    return env;
}

animation* make_animation(const std::string& name,
    const std::vector<float>& times, const std::vector<vec3f>& translation,
    const std::vector<vec4f>& rotation, const std::vector<vec3f>& scale,
    const std::vector<node*>& targets, bool bezier) {
    auto anm = new animation();
    anm->name = name;
    anm->times = times;
    anm->translation = translation;
    anm->rotation = rotation;
    anm->scale = scale;
    anm->type = (bezier) ? animation_type::bezier : animation_type::linear;
    anm->targets = targets;
    return anm;
}

scene* make_scene(const std::string& name, const std::vector<camera*>& cams,
    const std::vector<instance*>& ists, const std::vector<environment*>& envs) {
    auto scn = new scene();
    scn->name = name;
    scn->cameras = cams;
    scn->instances = ists;
    scn->environments = envs;
    auto add_elem = [](auto& elems, auto elem) {
        if (!elem) return;
        for (auto e : elems)
            if (e == elem) return;
        elems.push_back(elem);
    };
    for (auto ist : ists) {
        add_elem(scn->shapes, ist->shp);
        add_elem(scn->materials, ist->mat);
    }
    for (auto mat : scn->materials) {
        add_elem(scn->textures, mat->ke_txt.txt);
        add_elem(scn->textures, mat->kd_txt.txt);
        add_elem(scn->textures, mat->ks_txt.txt);
        add_elem(scn->textures, mat->norm_txt.txt);
    }
    for (auto env : scn->environments) {
        add_elem(scn->textures, env->ke_txt.txt);
    }
    return scn;
}

scene* make_scene(const std::string& name, const std::vector<camera*>& cams,
    const std::vector<instance*>& ists, const std::vector<environment*>& envs,
    const std::vector<node*>& ndes, const std::vector<animation*>& anms) {
    auto scn = make_scene(name, cams, ists, envs);
    scn->nodes = ndes;
    scn->animations = anms;
    return scn;
}

// example shapes
shape* make_floor_shape(const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_quad(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << tesselation, 1 << tesselation}, {size, size},
        {size / 2, size / 2});
    for (auto& p : shp->pos) p = {p.x, p.z, p.y};
    for (auto& n : shp->norm) n = {n.x, n.z, n.y};
    return shp;
}
shape* make_quad_shape(const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_quad(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << tesselation, 1 << tesselation}, {size, size},
        {size / 2, size / 2});
    return shp;
}
shape* make_sphere_shape(const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_sphere(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << (tesselation + 1), 1 << tesselation}, size, {size, size / 2});
    for (auto& p : shp->pos) p = {p.x, p.z, p.y};
    for (auto& n : shp->norm) n = {n.x, n.z, n.y};
    return shp;
}
shape* make_spherecube_shape(
    const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_sphere_cube(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        1 << tesselation, size, size / 2);
    return shp;
}
shape* make_sphereflipcap_shape(
    const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_sphere_flipcap(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << (tesselation + 1), 1 << tesselation}, size, {size, size / 2},
        {-0.75f, 0.75f});
    return shp;
}
shape* make_cube_shape(const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_cube(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << tesselation, 1 << tesselation, 1 << tesselation},
        {size, size, size}, {size / 2, size / 2, size / 2});
    return shp;
}
shape* make_cuberounded_shape(
    const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_cube_rounded(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << tesselation, 1 << tesselation, 1 << tesselation},
        {size, size, size}, {size / 2, size / 2, size / 2}, 0.15f * size);
    return shp;
}
shape* make_matball_shape(
    const std::string& name, int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_sphere(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << (tesselation + 1), 1 << tesselation}, size, {size, size / 2});
    for (auto& p : shp->pos) p = {p.x, p.z, p.y};
    for (auto& n : shp->norm) n = {n.x, n.z, n.y};
    return shp;
}
shape* make_quadstack_shape(const std::string& name, int stack_tesselation,
    int tesselation, float size) {
    auto shp = new shape();
    shp->name = name;
    make_quad_stack(shp->triangles, shp->pos, shp->norm, shp->texcoord,
        {1 << tesselation, 1 << tesselation, 1 << stack_tesselation},
        {size, size, size}, {size / 2, size / 2});
    return shp;
}
shape* make_cube_subdiv_shape(
    const std::string& name, int tesselation, float size) {
    auto pos = std::vector<vec3f>();
    auto norm = std::vector<vec3f>();
    auto quads = std::vector<vec4i>();
    make_cube(quads, pos, 0, size);
    for (auto i = 0; i < tesselation; i++) {
        subdivide_catmullclark(quads, pos);
    }
    if (tesselation) compute_normals(quads, pos, norm);
    auto shp = new shape();
    shp->name = name;
    shp->pos = pos;
    shp->norm = norm;
    shp->triangles = convert_quads_to_triangles(quads);
    return shp;
}
shape* make_suzanne_subdiv_shape(const std::string& name, int tesselation) {
    auto pos = std::vector<vec3f>();
    auto norm = std::vector<vec3f>();
    auto quads = std::vector<vec4i>();
    make_suzanne(quads, pos, 0);
    for (auto i = 0; i < tesselation; i++) {
        subdivide_catmullclark(quads, pos);
    }
    if (tesselation) compute_normals(quads, pos, norm);
    auto shp = new shape();
    shp->name = name;
    shp->pos = pos;
    shp->norm = norm;
    shp->triangles = convert_quads_to_triangles(quads);
    return shp;
}
shape* make_fvcube_subdiv_shape(
    const std::string& name, int tesselation, float size) {
    std::vector<vec4i> quads_pos, quads_norm, quads_texcoord;
    std::vector<vec3f> pos, norm;
    std::vector<vec2f> texcoord;
    make_fvcube(quads_pos, pos, quads_norm, norm, quads_texcoord, texcoord, 0,
        size, size / 2);
    for (auto i = 0; i < tesselation; i++) {
        subdivide_catmullclark(quads_pos, pos);
        subdivide_catmullclark(quads_texcoord, texcoord);
    }
    quads_norm = quads_pos;
    if (tesselation) compute_normals(quads_pos, pos, norm);
    auto shp = new shape();
    shp->name = name;
    auto quads = std::vector<vec4i>();
    convert_face_varying(quads, shp->pos, shp->norm, shp->texcoord, quads_pos,
        quads_norm, quads_texcoord, pos, norm, texcoord);
    shp->triangles = convert_quads_to_triangles(quads);
    return shp;
}
shape* make_hairball_shape(const std::string& name, int hair_tesselation,
    int tesselation, float size, const vec2f& len, const vec2f& noise,
    const vec2f& clump, const vec2f& radius) {
    auto shp1 = new shape();
    make_sphere_cube(shp1->triangles, shp1->pos, shp1->norm, shp1->texcoord,
        1 << 4, size * 0.8f, 1);
    auto shp = new shape();
    shp->name = name;
    make_hair(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius,
        {1 << tesselation, 1 << hair_tesselation}, shp1->triangles, shp1->pos,
        shp1->norm, shp1->texcoord, len, radius, noise, clump);
    delete shp1;
    return shp;
}

// example materials
material* make_emission_material(
    const std::string& name, const vec3f& col, texture* txt, texture* norm) {
    auto mat = make_material(name, zero3f, zero3f, 1);
    mat->ke = col;
    mat->ke_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}
material* make_matte_material(
    const std::string& name, const vec3f& col, texture* txt, texture* norm) {
    auto mat = make_material(name, col, zero3f, 1);
    mat->kd_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}
material* make_plastic_material(const std::string& name, const vec3f& col,
    float rs, texture* txt, texture* norm) {
    auto mat = make_material(name, col, {0.04f, 0.04f, 0.04f}, rs);
    mat->kd_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}
material* make_metal_material(const std::string& name, const vec3f& col,
    float rs, texture* txt, texture* norm) {
    auto mat = make_material(name, {0, 0, 0}, col, rs);
    mat->ks_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}
material* make_glass_material(const std::string& name, const vec3f& col,
    float rs, texture* txt, texture* norm) {
    auto mat = make_material(name, zero3f, {0.04f, 0.04f, 0.04f}, rs);
    mat->kt = col;
    mat->kt_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}
material* make_transparent_material(const std::string& name, const vec3f& col,
    float op, texture* txt, texture* norm) {
    auto mat = make_material(name, col, zero3f, 1);
    mat->op = op;
    mat->kd_txt.txt = txt;
    mat->norm_txt.txt = norm;
    return mat;
}

// example textures
texture* make_uvgrid_texture(const std::string& name, int res, int tile) {
    return make_texture(
        name, name + ".png", res, res, make_uvgrid_image(res, res, tile));
}
texture* make_grid_texture(const std::string& name, int res, int tile) {
    return make_texture(
        name, name + ".png", res, res, make_grid_image(res, res, tile));
}
texture* make_bump_texture(const std::string& name, int res, int tile) {
    return make_texture(
        name, name + ".png", res, res, make_bumpdimple_image(res, res, tile));
}
texture* make_bumpnorm_texture(
    const std::string& name, int res, int tile, float scale) {
    return make_texture(name, name + ".png", res, res,
        bump_to_normal_map(
            res, res, make_bumpdimple_image(res, res, tile), scale));
}
texture* make_sky_texture(const std::string& name, int res, float skyangle) {
    return make_texture(name, name + ".hdr", res * 2, res, {},
        make_sunsky_image(res * 2, res, skyangle));
}
texture* make_lights_texture(const std::string& name, int res, const vec3f& le,
    int nlights, float langle, float lwidth, float lheight) {
    return make_texture(name, name + ".hdr", res * 2, res, {},
        make_lights_image(res * 2, res, le, nlights, langle, lwidth, lheight));
}

// makes the cornell box scene
// http://graphics.cs.williams.edu/data
// http://www.graphics.cornell.edu/online/box/data.html
scene* make_cornellbox_scene(const std::string& name, bool envlight) {
    auto add_quad = [](scene* scn, std::string name, material* mat, vec3f pos,
                        vec3f rot = {0, 0, 0}, vec2f size = {2, 2}) {
        auto shp = new shape();
        shp->name = name;
        ygl::make_quad(shp->triangles, shp->pos, shp->norm, shp->texcoord,
            {1, 1}, size, {1, 1});
        scn->shapes.push_back(shp);
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->mat = mat;
        ist->frame = translation_frame(pos);
        if (rot != zero3f) {
            ist->frame = ist->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot.z * pi / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot.y * pi / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot.x * pi / 180);
        }
        scn->instances.push_back(ist);
    };

    auto add_box = [](scene* scn, std::string name, material* mat, vec3f pos,
                       vec3f rot = {0, 0, 0}, vec3f size = {2, 2, 2}) {
        auto shp = new shape();
        shp->name = name;
        shp->name = name;
        make_cube(shp->triangles, shp->pos, shp->norm, shp->texcoord, {1, 1, 1},
            size, {1, 1, 1});
        scn->shapes.push_back(shp);
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->mat = mat;
        ist->frame = translation_frame(pos);
        if (rot != zero3f) {
            ist->frame = ist->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot.z * pi / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot.y * pi / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot.x * pi / 180);
        }
        scn->instances.push_back(ist);
    };

    auto scn = new scene();
    scn->name = name;
    scn->cameras.push_back(
        make_camera("cam", {0, 1, 5.15f}, {0, 1, 0}, 27 * pi / 180, 1));
    scn->materials.push_back(make_material("white", {0.725f, 0.71f, 0.68f}));
    scn->materials.push_back(make_material("red", {0.63f, 0.065f, 0.05f}));
    scn->materials.push_back(make_material("green", {0.14f, 0.45f, 0.091f}));
    add_quad(scn, "floor", scn->materials[0], {0, 0, 0}, {-90, 0, 0});
    add_quad(scn, "ceiling", scn->materials[0], {0, 2, 0}, {90, 0, 0});
    add_quad(scn, "back", scn->materials[0], {0, 1, -1}, {0, 0, 0});
    add_quad(scn, "left", scn->materials[2], {+1, 1, 0}, {0, -90, 0});
    add_quad(scn, "right", scn->materials[1], {-1, 1, 0}, {0, 90, 0});
    add_box(scn, "tallbox", scn->materials[0], {-0.33f, 0.6f, -0.29f},
        {0, 15, 0}, {0.6f, 1.2f, 0.6f});
    add_box(scn, "shortbox", scn->materials[0], {0.33f, 0.3f, 0.33f},
        {0, -15, 0}, {0.6f, 0.6f, 0.6f});
    if (!envlight) {
        scn->materials.push_back(make_material("light", zero3f));
        scn->materials.back()->ke = {17, 12, 4};
        add_quad(scn, "light", scn->materials[3], {0, 1.999f, 0}, {90, 0, 0},
            {0.5f, 0.5f});
    } else {
        scn->environments.push_back(make_environment("env"));
    }
    return scn;
}

// make simple floor instances
instance* make_simple_floor(bool notexture = false) {
    return make_instance("floor", make_floor_shape("floor"),
        make_matte_material("floor",
            (notexture) ? vec3f{0.2f, 0.2f, 0.2f} : vec3f{1, 1, 1},
            (notexture) ? nullptr : make_grid_texture("grid")));
}

std::vector<instance*> make_simple_arealights() {
    return std::vector<instance*>{
        make_instance("light1", make_quad_shape("light1", 0, 4),
            make_emission_material("light1", {20, 20, 20}),
            lookat_frame({-4, 8, 8}, {0, 1, 0}, {0, 1, 0}, true)),
        make_instance("light2", make_quad_shape("light2", 0, 4),
            make_emission_material("light2", {20, 20, 20}),
            lookat_frame({+4, 8, 8}, {0, 1, 0}, {0, 1, 0}, true)),
    };
}

scene* make_simple_scene(const std::string& name,
    const std::vector<shape*>& shps, const std::vector<material*>& mats,
    bool envlight, const std::vector<animation*>& anms) {
    auto frames = std::vector<frame3f>();
    if (shps.size() >= 3) {
        frames = std::vector<frame3f>{
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-2.50f, 1, 0}},
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}},
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {+2.50f, 1, 0}}};
    } else if (shps.size() == 2) {
        frames = std::vector<frame3f>{
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1.25f, 1, 0}},
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {+1.25f, 1, 0}}};
    } else if (shps.size() == 1) {
        frames =
            std::vector<frame3f>{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}}};
    } else {
        throw std::runtime_error("number of shapes not supported");
    }

    auto cams = std::vector<camera*>();
    auto ists = std::vector<instance*>();
    auto envs = std::vector<environment*>();
    auto objs = std::vector<instance*>();

    if (shps.size() >= 3) {
        cams.push_back(make_camera(
            "cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 2.35f / 1.0f));
    } else if (shps.size() == 2) {
        cams.push_back(make_camera(
            "cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 16.0f / 9.0f));
    } else if (shps.size() == 1) {
        cams.push_back(
            make_camera("cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 1.0f));
    }

    ists.push_back(make_simple_floor());
    for (auto i = 0; i < shps.size(); i++) {
        ists.push_back(
            make_instance(shps[i]->name, shps[i], mats[i], frames[i % 3]));
        objs.push_back(ists.back());
    }

    if (envlight) {
        envs.push_back(
            make_environment("env", {1, 1, 1}, make_sky_texture("sky")));
    } else {
        for (auto lgt : make_simple_arealights()) ists.push_back(lgt);
    }

    if (anms.empty()) return make_scene(name, cams, ists, envs);

    auto ndes = std::vector<node*>();
    if (!anms.empty()) {
        for (auto cam : cams) {
            ndes.push_back(
                make_node(cam->name, cam, nullptr, nullptr, cam->frame));
        }
        for (auto ist : ists) {
            ndes.push_back(
                make_node(ist->name, nullptr, ist, nullptr, ist->frame));
            for (auto i = 0; i < anms.size(); i++) {
                if (objs[i] == ist) {
                    anms[i]->targets.push_back(ndes.back());
                    ndes.back()->frame = identity_frame3f;
                    ndes.back()->translation = frames[i % 3].o;
                }
            }
        }
        for (auto env : envs) {
            ndes.push_back(
                make_node(env->name, nullptr, nullptr, env, env->frame));
        }
    }

    return make_scene(name, cams, ists, envs, ndes, anms);
}

// Make a simple scene with single object.
scene* make_simple_scene(
    const std::string& name, shape* shp, material* mat, bool envlight) {
    auto cams = std::vector<camera*>();
    auto ists = std::vector<instance*>();
    auto envs = std::vector<environment*>();
    cams.push_back(
        make_camera("cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 1.0f));
    ists.push_back(make_simple_floor());
    ists.push_back(
        make_instance(shp->name, shp, mat, translation_frame({0, 1, 0})));
    if (envlight) {
        envs.push_back(
            make_environment("env", {1, 1, 1}, make_sky_texture("env")));
    } else {
        for (auto lgt : make_simple_arealights()) ists.push_back(lgt);
    }
    return make_scene(name, cams, ists, envs);
}

// Make a simple scene with single object.
scene* make_shape_scene(
    const std::string& name, shape* shp, material* mat, bool envlight) {
    auto cams = std::vector<camera*>();
    cams.push_back(
        make_camera("cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 1.0f));
    auto ists = std::vector<instance*>();
    auto envs = std::vector<environment*>();
    ists.push_back(
        make_instance("obj", shp, mat, translation_frame({0, 1, 0})));
    if (envlight) {
        envs.push_back(
            make_environment("env", {1, 1, 1}, make_sky_texture("sky")));
    } else {
        for (auto lgt : make_simple_arealights()) ists.push_back(lgt);
    }
    return make_scene(name, cams, ists, envs);
}

// Make a simple scene with single environment and a mirror ball.
scene* make_environment_scene(const std::string& name, environment* env) {
    auto cams = std::vector<camera*>();
    cams.push_back(
        make_camera("cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 1.0f));
    auto ists = std::vector<instance*>();
    auto envs = std::vector<environment*>();
    ists.push_back(make_instance("obj", make_sphere_shape("obj"),
        make_metal_material("obj", {1, 1, 1}, 0), identity_frame3f));
    envs.push_back(env);
    return make_scene(name, cams, ists, envs);
}

// instances shared functions
scene* make_random_instances_scene(const std::string& name, const vec2i& num,
    const bbox3f& bbox, uint64_t seed) {
    auto rscale = 0.9f * 0.25f *
                  min((bbox.max.x - bbox.min.x) / num.x,
                      (bbox.max.x - bbox.min.x) / num.y);

    auto cam = make_camera(
        "cam", {0, 5, 14}, {0, 1, 0}, 14.0f * pi / 180, 2.35f / 1.0f);

    auto ists = std::vector<instance*>();
    ists.push_back(make_instance("floor", make_floor_shape("floor"),
        make_matte_material("floor", {0.2f, 0.2f, 0.2f})));

    auto shps = std::vector<shape*>();
    auto mats = std::vector<material*>();

    auto cols = std::vector<vec3f>{
        {0.5f, 0.2f, 0.2f}, {0.2f, 0.5f, 0.2f}, {0.2f, 0.2f, 0.5f}};
    auto sid = 0;
    for (auto shp : {0, 1, 2}) {
        for (auto col : cols) {
            auto name = "shp" + std::to_string(sid++);
            mats.push_back(make_material(name, col));
            switch (shp) {
                case 0:
                    shps.push_back(make_sphere_shape(name, 4, 2 * rscale));
                    break;
                case 1:
                    shps.push_back(
                        make_sphereflipcap_shape(name, 4, 2 * rscale));
                    break;
                case 2:
                    shps.push_back(make_cuberounded_shape(name, 4, 2 * rscale));
                    break;
            }
        }
    }

    auto rng = make_rng(seed, 7);
    for (auto j = 0; j < num.y; j++) {
        for (auto i = 0; i < num.x; i++) {
            auto name = "ist" + std::to_string(j * num.x + i + 1);
            auto rpos = rand2f(rng);
            auto pos = vec3f{
                bbox.min.x + (bbox.max.x - bbox.min.x) *
                                 (i + 0.45f + 0.1f * rpos.x) / num.x,
                rscale,
                bbox.min.y + (bbox.max.y - bbox.min.y) *
                                 (j + 0.45f + 0.1f * rpos.y) / num.y,
            };
            auto idx = rand1i(rng, (int)mats.size());
            ists.push_back(make_instance(
                name, shps[idx], mats[idx], translation_frame(pos)));
        }
    }

    auto lpos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
    for (auto i = 0; i < 2; i++) {
        auto name = "light" + std::to_string(i + 1);
        ists.push_back(make_instance(name, make_quad_shape(name),
            make_emission_material(name, {80, 80, 80}),
            translation_frame(lpos[i])));
    }

    return make_scene(name, {cam}, ists, {});
}

}  // namespace ygl
