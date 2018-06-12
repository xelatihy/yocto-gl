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

#include <array>
#include <unordered_map>

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute shape normals. Supports only non-facevarying shapes.
void update_normals(const std::shared_ptr<shape>& shp) {
    if (!shp->triangles.empty()) {
        compute_normals(shp->triangles, shp->pos, shp->norm);
    } else if (!shp->lines.empty()) {
        compute_tangents(shp->lines, shp->pos, shp->norm);
    } else {
        shp->norm.assign(shp->pos.size(), {0, 0, 1});
    }
}

// Computes a shape bounding box.
void update_bbox(const std::shared_ptr<shape>& shp) {
    shp->bbox = invalid_bbox3f;
    for (auto p : shp->pos) shp->bbox += p;
}

// Updates the scene and scene's instances bounding boxes
void update_bbox(const std::shared_ptr<scene>& scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bbox(shp);
    }
    scn->bbox = invalid_bbox3f;
    for (auto ist : scn->instances) {
        ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);
        scn->bbox += ist->bbox;
    }
}

// Updates tesselation.
void update_tesselation(
    const std::shared_ptr<subdiv>& sbd, std::shared_ptr<shape> shp) {
    shp->name = sbd->name;
    auto quads_pos = sbd->quads_pos;
    auto quads_texcoord = sbd->quads_texcoord;
    auto quads_color = sbd->quads_color;
    auto pos = sbd->pos;
    auto texcoord = sbd->texcoord;
    auto color = sbd->color;
    subdivide_catmullclark(quads_pos, pos, sbd->level);
    subdivide_catmullclark(quads_texcoord, texcoord, sbd->level, true);
    subdivide_catmullclark(quads_color, color, sbd->level);
    auto norm = std::vector<vec3f>();
    if (sbd->compute_normals) compute_normals(quads_pos, pos, norm);
    auto quads = quads_pos;
    convert_face_varying(quads, shp->pos, shp->norm, shp->texcoord, shp->color,
        quads_pos, quads_pos, quads_texcoord, quads_color, pos, norm, texcoord,
        color);
    shp->triangles = convert_quads_to_triangles(quads);
    update_bbox(shp);
}
void update_tesselation(const std::shared_ptr<scene>& scn) {
    for (auto ist : scn->instances) {
        if (!ist->sbd) continue;
        update_tesselation(ist->sbd, ist->shp);
    }
}

// Update animation transforms
void update_transforms(const std::shared_ptr<animation>& anm, float time,
    const std::string& anim_group) {
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
void update_transforms(const std::shared_ptr<node>& nde,
    const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame(nde->translation) *
                 rotation_frame(nde->rotation) * scaling_frame(nde->scale);
    if (nde->ist) nde->ist->frame = frame;
    if (nde->cam) nde->cam->frame = frame;
    if (nde->env) nde->env->frame = frame;
    for (auto& child : nde->children) update_transforms(child.lock(), frame);
}

// Update node transforms
void update_transforms(const std::shared_ptr<scene>& scn, float time,
    const std::string& anim_group) {
    for (auto& agr : scn->animations) update_transforms(agr, time, anim_group);
    for (auto& nde : scn->nodes) nde->children.clear();
    for (auto& nde : scn->nodes)
        if (nde->parent) nde->parent->children.push_back(nde);
    for (auto& nde : scn->nodes)
        if (!nde->parent) update_transforms(nde);
}

// Compute animation range
vec2f compute_animation_range(
    const std::shared_ptr<scene>& scn, const std::string& anim_group) {
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
void update_lights(
    const std::shared_ptr<scene>& scn, bool do_shapes, bool do_environments) {
    if (do_shapes) {
        for (auto shp : scn->shapes) shp->elem_cdf.clear();
    }
    if (do_environments) {
        for (auto env : scn->environments) env->elem_cdf.clear();
    }
    scn->lights.clear();

    for (auto ist : scn->instances) {
        if (!ist->mat || ist->mat->ke == zero3f) continue;
        if (ist->shp->triangles.empty()) continue;
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
void update_shape_cdf(const std::shared_ptr<shape>& shp) {
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
void update_environment_cdf(std::shared_ptr<environment> env) {
    env->elem_cdf.clear();
    auto txt = env->ke_txt;
    if (!txt) return;
    env->elem_cdf.resize(txt->img.width() * txt->img.height());
    if (!txt->img.empty()) {
        for (auto i = 0; i < env->elem_cdf.size(); i++) {
            auto th = (i / txt->img.width() + 0.5f) * pi / txt->img.height();
            env->elem_cdf[i] = max(xyz(txt->img[i])) * sin(th);
            if (i) env->elem_cdf[i] += env->elem_cdf[i - 1];
        }
    } else {
        throw std::runtime_error("empty texture");
    }
}

// Build a shape BVH
void update_bvh(const std::shared_ptr<shape>& shp, bool equalsize) {
    if (!shp->bvh) shp->bvh = std::make_shared<bvh_tree>();
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    shp->bvh->points = shp->points;
    shp->bvh->lines = shp->lines;
    shp->bvh->triangles = shp->triangles;
    build_bvh(shp->bvh, equalsize);
}

// Build a scene BVH
void update_bvh(
    const std::shared_ptr<scene>& scn, bool do_shapes, bool equalsize) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bvh(shp, equalsize);
    }

    // tree bvh
    if (!scn->bvh) scn->bvh = std::make_shared<bvh_tree>();
    scn->bvh->ist_frames.resize(scn->instances.size());
    scn->bvh->ist_inv_frames.resize(scn->instances.size());
    scn->bvh->ist_bvhs.resize(scn->instances.size());
    for (auto i = 0; i < scn->instances.size(); i++) {
        auto ist = scn->instances[i];
        scn->bvh->ist_frames[i] = ist->frame;
        scn->bvh->ist_inv_frames[i] = inverse(ist->frame, false);
        scn->bvh->ist_bvhs[i] = ist->shp->bvh;
    }
    build_bvh(scn->bvh, equalsize);
}

// Refits a scene BVH
void refit_bvh(const std::shared_ptr<shape>& shp) {
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    refit_bvh(shp->bvh);
}

// Refits a scene BVH
void refit_bvh(const std::shared_ptr<scene>& scn, bool do_shapes) {
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
void add_missing_names(const std::shared_ptr<scene>& scn) {
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
void add_missing_normals(const std::shared_ptr<scene>& scn) {
    for (auto shp : scn->shapes) {
        if (shp->norm.empty()) update_normals(shp);
    }
}

// Add missing tangent space if needed.
void add_missing_tangent_space(const std::shared_ptr<scene>& scn) {
    for (auto ist : scn->instances) {
        if (!ist->shp->tangsp.empty() || ist->shp->texcoord.empty()) continue;
        if (!ist->mat || (!ist->mat->norm_txt && !ist->mat->bump_txt)) continue;
        if (!ist->shp->triangles.empty()) {
            if (ist->shp->norm.empty()) update_normals(ist->shp);
            compute_tangent_space(ist->shp->triangles, ist->shp->pos,
                ist->shp->norm, ist->shp->texcoord, ist->shp->tangsp);
        } else {
            throw std::runtime_error("type not supported");
        }
    }
}

// Checks for validity of the scene.
std::vector<std::string> validate(
    const std::shared_ptr<scene>& scn, bool skip_textures) {
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
    auto check_empty_textures =
        [&errs](const std::vector<std::shared_ptr<texture>>& vals) {
            for (auto val : vals) {
                if (val->img.empty())
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
    const std::shared_ptr<scene>& scn, const ray3f& ray, bool find_any) {
    auto iid = 0;
    auto isec = scene_intersection();
    if (!intersect_bvh(
            scn->bvh, ray, find_any, isec.dist, iid, isec.ei, isec.uv))
        return {};
    isec.ist = scn->instances[iid];
    return isec;
}

// Shape element normal.
vec3f eval_elem_norm(const std::shared_ptr<shape>& shp, int ei) {
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

// Shape element normal.
vec4f eval_elem_tangsp(const std::shared_ptr<shape>& shp, int ei) {
    auto tangsp = zero4f;
    if (!shp->triangles.empty()) {
        auto t = shp->triangles[ei];
        auto norm =
            triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        auto txty = std::pair<vec3f, vec3f>();
        if (shp->texcoord.empty()) {
            txty = triangle_tangents_fromuv(shp->pos[t.x], shp->pos[t.y],
                shp->pos[t.z], {0, 0}, {1, 0}, {0, 1});
        } else {
            txty = triangle_tangents_fromuv(shp->pos[t.x], shp->pos[t.y],
                shp->pos[t.z], shp->texcoord[t.x], shp->texcoord[t.y],
                shp->texcoord[t.z]);
        }
        auto tx = txty.first, ty = txty.second;
        tx = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        tangsp = {tx.x, tx.y, tx.z, s};
    }
    return tangsp;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T eval_elem(const std::shared_ptr<shape>& shp, const std::vector<T>& vals,
    int ei, const vec2f& uv) {
    if (vals.empty()) return {};
    if (!shp->triangles.empty()) {
        return interpolate_triangle(vals, shp->triangles[ei], uv);
    } else if (!shp->lines.empty()) {
        return interpolate_line(vals, shp->lines[ei], uv.x);
    } else if (!shp->points.empty()) {
        return vals[shp->points[ei]];
    } else if (!shp->pos.empty()) {
        return vals[ei];
    } else {
        return {};
    }
}

// Shape values interpolated using barycentric coordinates
vec3f eval_pos(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->pos, ei, uv);
}
vec3f eval_norm(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->norm.empty()) return eval_elem_norm(shp, ei);
    return normalize(eval_elem(shp, shp->norm, ei, uv));
}
vec2f eval_texcoord(
    const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->texcoord.empty()) return uv;
    return eval_elem(shp, shp->texcoord, ei, uv);
}
vec4f eval_color(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->color.empty()) return {1, 1, 1, 1};
    return eval_elem(shp, shp->color, ei, uv);
}
float eval_radius(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->radius.empty()) return 0.001f;
    return eval_elem(shp, shp->radius, ei, uv);
}
vec4f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->tangsp.empty()) return eval_elem_tangsp(shp, ei);
    return eval_elem(shp, shp->tangsp, ei, uv);
}
vec3f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv,
    bool& left_handed) {
    auto tangsp = (shp->tangsp.empty()) ? eval_elem_tangsp(shp, ei) :
                                          eval_elem(shp, shp->tangsp, ei, uv);
    left_handed = tangsp.w < 0;
    return {tangsp.x, tangsp.y, tangsp.z};
}

// Instance values interpolated using barycentric coordinates.
vec3f eval_pos(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return transform_point(ist->frame, eval_pos(ist->shp, ei, uv));
}
vec3f eval_norm(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return transform_direction(ist->frame, eval_norm(ist->shp, ei, uv));
}
vec2f eval_texcoord(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_texcoord(ist->shp, ei, uv);
}
vec4f eval_color(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_color(ist->shp, ei, uv);
}
float eval_radius(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_radius(ist->shp, ei, uv);
}
vec3f eval_tangsp(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv,
    bool& left_handed) {
    return transform_direction(
        ist->frame, eval_tangsp(ist->shp, ei, uv, left_handed));
}
// Instance element values.
vec3f eval_elem_norm(const std::shared_ptr<instance>& ist, int ei) {
    return transform_direction(ist->frame, eval_elem_norm(ist->shp, ei));
}
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv, vec3f o) {
    if (!ist->shp->triangles.empty()) {
        auto n = eval_norm(ist, ei, uv);
        if (ist->mat && ist->mat->norm_txt) {
            auto texcoord = eval_texcoord(ist, ei, uv);
            auto left_handed = false;
            auto txt = xyz(eval_texture(ist->mat->norm_txt, texcoord));
            txt = txt * 2 - vec3f{1, 1, 1};
            txt.y = -txt.y;  // flip vertical axis to align green with image up
            auto tu = orthonormalize(eval_tangsp(ist, ei, uv, left_handed), n);
            auto tv = normalize(cross(n, tu) * (left_handed ? -1.0f : 1.0f));
            n = normalize(txt.x * tu + txt.y * tv + txt.z * n);
        }
        if (ist->mat && ist->mat->double_sided && dot(n, o) < 0) n = -n;
        return n;
    } else if (!ist->shp->lines.empty()) {
        return orthonormalize(o, eval_norm(ist, ei, uv));
    } else {
        return o;
    }
}

// Environment texture coordinates from the direction.
vec2f eval_texcoord(const std::shared_ptr<environment>& env, vec3f w) {
    auto wl = transform_direction_inverse(env->frame, w);
    auto uv = vec2f{
        atan2(wl.z, wl.x) / (2 * pi), acos(clamp(wl.y, -1.0f, 1.0f)) / pi};
    if (uv.x < 0) uv.x += 1;
    return uv;
}
// Evaluate the environment direction.
vec3f eval_direction(const std::shared_ptr<environment>& env, const vec2f& uv) {
    return transform_direction(
        env->frame, {cos(uv.x * 2 * pi) * sin(uv.y * pi), cos(uv.y * pi),
                        sin(uv.x * 2 * pi) * sin(uv.y * pi)});
}
// Evaluate the environment color.
vec3f eval_environment(const std::shared_ptr<environment>& env, vec3f w) {
    auto ke = env->ke;
    if (env->ke_txt) {
        ke *= xyz(eval_texture(env->ke_txt, eval_texcoord(env, w)));
    }
    return ke;
}

// Evaluate a texture
vec4f eval_texture(const std::shared_ptr<texture>& txt, const vec2f& texcoord) {
    if (!txt || txt->img.empty()) return {1, 1, 1, 1};

    // get image width/height
    auto w = txt->img.width(), h = txt->img.height();

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (txt->clamp) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * w;
        t = clamp(texcoord.y, 0.0f, 1.0f) * h;
    } else {
        s = std::fmod(texcoord.x, 1.0f) * w;
        if (s < 0) s += w;
        t = std::fmod(texcoord.y, 1.0f) * h;
        if (t < 0) t += h;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, w - 1), j = clamp((int)t, 0, h - 1);
    auto ii = (i + 1) % w, jj = (j + 1) % h;
    auto u = s - i, v = t - j;

    // handle interpolation
    return txt->img[{i, j}] * (1 - u) * (1 - v) +
           txt->img[{i, jj}] * (1 - u) * v + txt->img[{ii, j}] * u * (1 - v) +
           txt->img[{ii, jj}] * u * v;
}

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const std::shared_ptr<camera>& cam) {
    return 2 * std::atan(cam->height / (2 * cam->focal));
}
void set_camera_fovy(
    const std::shared_ptr<camera>& cam, float fovy, float aspect, float width) {
    cam->width = width;
    cam->height = width / aspect;
    cam->focal = cam->height / (2 * std::tan(fovy / 2));
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera_ray(
    const std::shared_ptr<camera>& cam, const vec2f& uv, const vec2f& luv) {
    auto dist = cam->focal;
    if (cam->focus < flt_max) {
        dist = cam->focal * cam->focus / (cam->focus - cam->focal);
    }
    auto e = vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    // auto q = vec3f{cam->width * (uv.x - 0.5f),
    //     cam->height * (uv.y - 0.5f), dist};
    // X flipped for mirror
    auto q =
        vec3f{cam->width * (0.5f - uv.x), cam->height * (uv.y - 0.5f), dist};
    auto ray = make_ray(transform_point(cam->frame, e),
        transform_direction(cam->frame, normalize(e - q)));
    return ray;
}

// Generates a ray from a camera for pixel coordinates `ij`, the
// resolution `res`, the sub-pixel coordinates `puv` and the lens
// coordinates `luv` and the image resolution `res`.
ray3f eval_camera_ray(const std::shared_ptr<camera>& cam, const vec2i& ij,
    const vec2i& imsize, const vec2f& puv, const vec2f& luv) {
    auto uv = vec2f{(ij.x + puv.x) / imsize.x, (ij.y + puv.y) / imsize.y};
    return eval_camera_ray(cam, uv, luv);
}

// Evaluates material parameters.
vec3f eval_emission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->ke * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->ke_txt, eval_texcoord(ist, ei, uv)));
}
vec3f eval_diffuse(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).z;
        return kb * (1 - km);
    }
}
vec3f eval_specular(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->ks * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).z;
        return kb * km + vec3f{0.04f, 0.04f, 0.04f} * (1 - km);
    }
}
float eval_roughness(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    if (!ist->mat->base_metallic) {
        if (!ist->mat->gltf_textures) {
            auto rs =
                ist->mat->rs *
                eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).x;
            return rs * rs;
        } else {
            auto gs =
                (1 - ist->mat->rs) *
                eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).w;
            auto rs = 1 - gs;
            return rs * rs;
        }
    } else {
        auto rs = ist->mat->rs *
                  eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).y;
        return rs * rs;
    }
}
vec3f eval_transmission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->kt * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->kt_txt, eval_texcoord(ist, ei, uv)));
}
float eval_opacity(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    return ist->mat->op * eval_color(ist->shp, ei, uv).w *
           eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)).w *
           eval_texture(ist->mat->op_txt, eval_texcoord(ist, ei, uv)).x;
}

// Evaluates the bsdf at a location.
bsdf eval_bsdf(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    auto f = bsdf();
    f.kd = eval_diffuse(ist, ei, uv);
    f.ks = eval_specular(ist, ei, uv);
    f.kt = eval_transmission(ist, ei, uv);
    f.rs = eval_roughness(ist, ei, uv);
    f.refract = (ist && ist->mat) ? ist->mat->refract : false;
    if (f.kd != zero3f) {
        f.rs = clamp(f.rs, 0.03f * 0.03f, 1.0f);
    } else if (f.rs <= 0.03f * 0.03f)
        f.rs = 0;
    return f;
}
bool is_delta_bsdf(const bsdf& f) { return f.rs == 0 && f.kd == zero3f; }

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const std::shared_ptr<shape>& shp, float re, const vec2f& ruv) {
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

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Merge scene into one another
void merge_into(const std::shared_ptr<scene>& merge_into,
    std::shared_ptr<scene> merge_from) {
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

void print_stats(const std::shared_ptr<scene>& scn) {
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

    uint64_t texel_hdr = 0;
    uint64_t texel_ldr = 0;

    uint64_t memory_imgs = 0;
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

    for (auto txt : scn->textures) { texel_hdr += txt->img.size(); }
    memory_imgs = texel_hdr * sizeof(vec4f) + texel_ldr * sizeof(vec4b);

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
    println("texel_hdr: {}", texel_hdr);
    println("texel_ldr: {}", texel_ldr);
    println("memory_imgs: {}", memory_imgs);
    println("memory_elems: {}", memory_elems);
    println("memory_verts: {}", memory_verts);
    println("bbox_scn: {} {}", bbox.min, bbox.max);
    println("bbox_min   : {}", bbox.min);
    println("bbox_max   : {}", bbox.max);
    println("bbox_size  : {}", bbox.max - bbox.min);
    println("bbox_center: {}", (bbox.max + bbox.min) / 2);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

std::shared_ptr<camera> make_camera(const std::string& name,
    const frame3f& frame, float width, float height, float focal, float focus,
    float aperture) {
    auto cam = std::make_shared<camera>();
    cam->name = name;
    cam->frame = frame;
    cam->width = width;
    cam->height = height;
    cam->focal = focal;
    cam->focus = focus;
    cam->aperture = aperture;
    cam->near = 0.01f;
    cam->far = flt_max;
    return cam;
};

// add missing camera
std::shared_ptr<camera> make_bbox_camera(const std::string& name,
    const bbox3f& bbox, float width, float height, float focal) {
    auto bbox_center = (bbox.max + bbox.min) / 2.0f;
    auto bbox_size = bbox.max - bbox.min;
    auto bbox_msize = max(bbox_size.x, max(bbox_size.y, bbox_size.z));
    auto cam = std::make_shared<camera>();
    cam->name = name;
    auto camera_dir = vec3f{1, 0.4f, 1};
    auto from = camera_dir * bbox_msize + bbox_center;
    auto to = bbox_center;
    auto up = vec3f{0, 1, 0};
    cam->frame = lookat_frame(from, to, up);
    cam->ortho = false;
    cam->width = width;
    cam->height = height;
    cam->focal = focal;
    cam->focus = flt_max;
    cam->aperture = 0;
    return cam;
}

std::shared_ptr<material> make_material(
    const std::string& name, vec3f kd, vec3f ks, float rs) {
    auto mat = std::make_shared<material>();
    mat->name = name;
    mat->kd = kd;
    mat->ks = ks;
    mat->rs = rs;
    return mat;
}

std::shared_ptr<texture> make_texture(
    const std::string& name, const std::string& path, const image4f& img) {
    auto txt = std::make_shared<texture>();
    txt->name = name;
    txt->path = path;
    txt->img = img;
    return txt;
}

std::shared_ptr<shape> make_shape(const std::string& name,
    const std::string& path, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius) {
    auto shp = std::make_shared<shape>();
    shp->name = name;
    shp->path = path;
    shp->pos = pos;
    shp->norm = norm;
    shp->texcoord = texcoord;
    shp->color = color;
    shp->radius = radius;
    shp->lines = lines;
    shp->triangles = triangles;
    return shp;
}

std::shared_ptr<subdiv> make_subdiv(const std::string& name,
    const std::string& path, int level, const std::vector<vec4i>& quads_pos,
    const std::vector<vec3f>& pos, const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec2f>& texcoord, const std::vector<vec4i>& quads_color,
    const std::vector<vec4f>& color) {
    auto sbd = std::make_shared<subdiv>();
    sbd->name = name;
    sbd->path = path;
    sbd->level = level;
    sbd->pos = pos;
    sbd->texcoord = texcoord;
    sbd->color = color;
    sbd->quads_pos = quads_pos;
    sbd->quads_texcoord = quads_texcoord;
    sbd->quads_color = quads_color;
    return sbd;
}

std::shared_ptr<instance> make_instance(const std::string& name,
    std::shared_ptr<shape> shp, std::shared_ptr<material> mat,
    std::shared_ptr<subdiv> sbd, const frame3f& frame) {
    auto ist = std::make_shared<instance>();
    ist->name = name;
    ist->shp = shp;
    ist->mat = mat;
    ist->sbd = sbd;
    ist->frame = frame;
    return ist;
}

std::shared_ptr<node> make_node(const std::string& name,
    std::shared_ptr<camera> cam, std::shared_ptr<instance> ist,
    std::shared_ptr<environment> env, const frame3f& frame) {
    auto nde = std::make_shared<node>();
    nde->name = name;
    nde->cam = cam;
    nde->ist = ist;
    nde->env = env;
    nde->frame = frame;
    return nde;
}

std::shared_ptr<environment> make_environment(const std::string& name, vec3f ke,
    std::shared_ptr<texture> ke_txt, const frame3f& frame) {
    auto env = std::make_shared<environment>();
    env->name = name;
    env->ke = ke;
    env->ke_txt = ke_txt;
    env->frame = frame;
    return env;
}

std::shared_ptr<animation> make_animation(const std::string& name,
    const std::string& path, const std::vector<float>& times,
    const std::vector<vec3f>& translation, const std::vector<vec4f>& rotation,
    const std::vector<vec3f>& scale,
    const std::vector<std::shared_ptr<node>>& targets, bool bezier) {
    auto anm = std::make_shared<animation>();
    anm->name = name;
    anm->path = path;
    anm->times = times;
    anm->translation = translation;
    anm->rotation = rotation;
    anm->scale = scale;
    anm->type = (bezier) ? animation_type::bezier : animation_type::linear;
    anm->targets = targets;
    return anm;
}

std::shared_ptr<scene> make_scene(const std::string& name,
    const std::vector<std::shared_ptr<camera>>& cams,
    const std::vector<std::shared_ptr<instance>>& ists,
    const std::vector<std::shared_ptr<environment>>& envs) {
    auto scn = std::make_shared<scene>();
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
        add_elem(scn->textures, mat->ke_txt);
        add_elem(scn->textures, mat->kd_txt);
        add_elem(scn->textures, mat->ks_txt);
        add_elem(scn->textures, mat->norm_txt);
    }
    for (auto env : scn->environments) { add_elem(scn->textures, env->ke_txt); }
    return scn;
}

std::shared_ptr<scene> make_scene(const std::string& name,
    const std::vector<std::shared_ptr<camera>>& cams,
    const std::vector<std::shared_ptr<instance>>& ists,
    const std::vector<std::shared_ptr<environment>>& envs,
    const std::vector<std::shared_ptr<node>>& ndes,
    const std::vector<std::shared_ptr<animation>>& anms) {
    auto scn = make_scene(name, cams, ists, envs);
    scn->nodes = ndes;
    scn->animations = anms;
    return scn;
}

}  // namespace ygl
