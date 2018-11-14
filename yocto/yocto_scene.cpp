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
#include "yocto_random.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes a shape bounding box.
bbox3f compute_shape_bounds(const yocto_shape& shape) {
    auto bbox = invalid_bbox3f;
    for (auto p : shape.positions) bbox += p;
    return bbox;
}

// Computes a surface bounding box.
bbox3f compute_surface_bounds(const yocto_surface& surface) {
    auto bbox = invalid_bbox3f;
    for (auto p : surface.positions) bbox += p;
    return bbox;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_scene_bounds(const yocto_scene& scene) {
    auto shape_bbox = vector<bbox3f>(scene.shapes.size());
    for (auto shape_id = 0; shape_id < scene.shapes.size(); shape_id++)
        shape_bbox[shape_id] = compute_shape_bounds(scene.shapes[shape_id]);
    auto surface_bbox = vector<bbox3f>(scene.shapes.size());
    for (auto surface_id = 0; surface_id < scene.surfaces.size(); surface_id++)
        surface_bbox[surface_id] = compute_surface_bounds(
            scene.surfaces[surface_id]);
    auto bbox = invalid_bbox3f;
    for (auto& instance : scene.instances) {
        if (instance.shape >= 0) {
            bbox += transform_bbox(instance.frame, shape_bbox[instance.shape]);
        } else if (instance.surface >= 0) {
            bbox += transform_bbox(
                instance.frame, surface_bbox[instance.surface]);
        }
    }
    return bbox;
}

// Compute vertex normals
vector<vec3f> compute_shape_normals(const yocto_shape& shape) {
    if (!shape.points.empty()) {
        return vector<vec3f>(shape.positions.size(), {0, 0, 1});
    } else if (!shape.lines.empty()) {
        return compute_vertex_tangents(shape.lines, shape.positions);
    } else if (!shape.triangles.empty()) {
        return compute_vertex_normals(shape.triangles, shape.positions);
    } else {
        return {};
    }
}

// Compute vertex normals
vector<vec3f> compute_surface_normals(const yocto_surface& shape) {
    if (!shape.quads_positions.empty()) {
        return compute_vertex_normals(shape.quads_positions, shape.positions);
    } else {
        return {};
    }
}

// Apply subdivision and displacement rules.
yocto_shape subdivide_shape(const yocto_shape& shape, int subdivision_level,
    bool catmull_clark, bool compute_normals) {
    if (!subdivision_level) return shape;
    auto subdivided = shape;
    if (!subdivided.points.empty()) {
        log_error("point subdivision not supported");
    } else if (!subdivided.lines.empty()) {
        for (auto l = 0; l < subdivision_level; l++) {
            auto lines                                  = subdivided.lines;
            tie(subdivided.lines, subdivided.positions) = subdivide_lines(
                lines, subdivided.positions);
            tie(ignore, subdivided.normals) = subdivide_lines(
                lines, subdivided.normals);
            tie(ignore, subdivided.texturecoords) = subdivide_lines(
                lines, subdivided.texturecoords);
            tie(ignore, subdivided.colors) = subdivide_lines(
                lines, subdivided.colors);
            tie(ignore, subdivided.radius) = subdivide_lines(
                lines, subdivided.radius);
        }
    } else if (!subdivided.triangles.empty()) {
        for (auto l = 0; l < subdivision_level; l++) {
            auto triangles = subdivided.triangles;
            tie(subdivided.triangles, subdivided.positions) = subdivide_triangles(
                triangles, subdivided.positions);
            tie(ignore, subdivided.normals) = subdivide_triangles(
                triangles, subdivided.normals);
            tie(ignore, subdivided.texturecoords) = subdivide_triangles(
                triangles, subdivided.texturecoords);
            tie(ignore, subdivided.colors) = subdivide_triangles(
                triangles, subdivided.colors);
        }
    } else if (!subdivided.quads.empty() && !catmull_clark) {
        for (auto l = 0; l < subdivision_level; l++) {
            auto quads                                  = subdivided.quads;
            tie(subdivided.quads, subdivided.positions) = subdivide_quads(
                quads, subdivided.positions);
            tie(ignore, subdivided.normals) = subdivide_quads(
                quads, subdivided.normals);
            tie(ignore, subdivided.texturecoords) = subdivide_quads(
                quads, subdivided.texturecoords);
            tie(ignore, subdivided.colors) = subdivide_quads(
                quads, subdivided.colors);
        }
    } else if (!subdivided.quads.empty() && catmull_clark) {
        for (auto l = 0; l < subdivision_level; l++) {
            auto quads                                  = subdivided.quads;
            tie(subdivided.quads, subdivided.positions) = subdivide_catmullclark(
                quads, subdivided.positions);
            tie(ignore, subdivided.normals) = subdivide_catmullclark(
                quads, subdivided.normals);
            tie(ignore, subdivided.texturecoords) = subdivide_catmullclark(
                quads, subdivided.texturecoords);
            tie(ignore, subdivided.colors) = subdivide_catmullclark(
                quads, subdivided.colors);
        }
    }

    if (compute_normals) {
        subdivided.normals = compute_shape_normals(subdivided);
    }

    return subdivided;
}
// Apply subdivision and displacement rules.
yocto_surface subdivide_surface(const yocto_surface& surface,
    int subdivision_level, bool catmull_clark, bool compute_normals) {
    if (!subdivision_level) return surface;
    auto subdivided    = surface;
    auto subdivide_ids = [](const vector<vec4i>& quads, const vector<int>& ids) {
        if (ids.empty()) return vector<int>{};
        auto new_ids = vector<int>();
        for (int quad_id = 0; quad_id < quads.size(); quad_id++) {
            auto quad = quads[quad_id];
            for (auto i = 0; i < (quad[2] == quad[3] ? 3 : 4); i++)
                new_ids.push_back(ids[quad_id]);
        }
        return new_ids;
    };
    if (!subdivided.quads_positions.empty() && !catmull_clark) {
        for (auto l = 0; l < subdivision_level; l++) {
            subdivided.quads_materials = subdivide_ids(
                subdivided.quads_positions, subdivided.quads_materials);
            tie(subdivided.quads_positions, subdivided.positions) = subdivide_quads(
                subdivided.quads_positions, subdivided.positions);
            tie(subdivided.quads_normals, subdivided.normals) = subdivide_quads(
                subdivided.quads_normals, subdivided.normals);
            tie(subdivided.quads_texturecoords, subdivided.texturecoords) = subdivide_quads(
                subdivided.quads_texturecoords, subdivided.texturecoords);
        }
    } else if (!subdivided.quads_positions.empty() && catmull_clark) {
        for (auto l = 0; l < subdivision_level; l++) {
            subdivided.quads_materials = subdivide_ids(
                subdivided.quads_positions, subdivided.quads_materials);
            tie(subdivided.quads_positions, subdivided.positions) = subdivide_catmullclark(
                subdivided.quads_positions, subdivided.positions);
            tie(subdivided.quads_texturecoords, subdivided.texturecoords) = subdivide_catmullclark(
                subdivided.quads_texturecoords, subdivided.texturecoords, true);
        }
    }

    if (compute_normals) {
        if (!subdivided.quads_positions.empty())
            subdivided.quads_normals = subdivided.quads_positions;
        subdivided.normals = compute_surface_normals(subdivided);
    }

    return subdivided;
}
yocto_shape displace_shape(const yocto_shape& shape,
    const yocto_texture& displacement, bool compute_normals) {
    if (shape.texturecoords.empty()) {
        log_error("missing texture coordinates");
        return shape;
    }
    auto displaced_shape = shape;
    auto normals = (shape.normals.empty()) ? compute_shape_normals(shape) :
                                             shape.normals;
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
        displaced_shape.positions[vid] += normals[vid] *
                                          displacement.height_scale *
                                          mean(make_shorter_vec(
                                              evaluate_texture(displacement,
                                                  shape.texturecoords[vid])));
    }

    if (compute_normals)
        displaced_shape.normals = compute_shape_normals(displaced_shape);

    return displaced_shape;
}
yocto_surface displace_surface(const yocto_surface& surface,
    const yocto_texture& displacement, bool compute_normals) {
    if (surface.texturecoords.empty()) {
        log_error("missing texture coordinates");
        return surface;
    }
    auto displaced_surface = surface;
    auto offset            = vector<float>(surface.positions.size(), 0);
    auto count             = vector<int>(surface.positions.size(), 0);
    for (auto fid = 0; fid < surface.quads_positions.size(); fid++) {
        auto qpos = surface.quads_positions[fid];
        auto qtxt = surface.quads_texturecoords[fid];
        for (auto i = 0; i < 4; i++) {
            offset[qpos[i]] += displacement.height_scale *
                               mean(make_shorter_vec(evaluate_texture(displacement,
                                   surface.texturecoords[qtxt[i]])));
            count[qpos[i]] += 1;
        }
    }
    auto normals = compute_vertex_normals(
        surface.quads_positions, surface.positions);
    for (auto vid = 0; vid < surface.positions.size(); vid++) {
        displaced_surface.positions[vid] += normals[vid] * offset[vid] /
                                            count[vid];
    }

    if (compute_normals) {
        displaced_surface.quads_normals = displaced_surface.quads_positions;
        displaced_surface.normals = compute_surface_normals(displaced_surface);
    }

    return displaced_surface;
}

// Updates tesselation.
void tesselate_shapes_and_surfaces(yocto_scene& scene) {
    auto scope = log_trace_scoped("tesselating surfaces");
    for (auto& shape : scene.shapes) {
        auto& material = scene.materials[shape.material];
        if (!shape.subdivision_level && material.displacement_texture < 0)
            continue;
        auto tesselated_shape = shape;
        if (shape.subdivision_level) {
            tesselated_shape = subdivide_shape(tesselated_shape,
                tesselated_shape.subdivision_level,
                tesselated_shape.catmull_clark,
                tesselated_shape.compute_vertex_normals);
            tesselated_shape.subdivision_level = 0;
        }
        if (material.displacement_texture >= 0) {
            tesselated_shape = displace_shape(tesselated_shape,
                scene.textures[material.displacement_texture],
                shape.compute_vertex_normals);
        }
        shape = tesselated_shape;
    }
    for (auto& surface : scene.surfaces) {
        auto& material = scene.materials[surface.materials.front()];
        if (!surface.subdivision_level && material.displacement_texture < 0)
            continue;
        auto tesselated_surface = surface;
        if (surface.subdivision_level) {
            tesselated_surface = subdivide_surface(tesselated_surface,
                tesselated_surface.subdivision_level,
                tesselated_surface.catmull_clark,
                tesselated_surface.compute_vertex_normals);
            tesselated_surface.subdivision_level = 0;
        }
        if (material.displacement_texture >= 0) {
            tesselated_surface = displace_surface(tesselated_surface,
                scene.textures[material.displacement_texture],
                surface.compute_vertex_normals);
        }
        surface = tesselated_surface;
    }
}

// Update animation transforms
void update_transforms(yocto_scene& scene, yocto_animation& animation,
    float time, const string& anim_group) {
    if (anim_group != "" && anim_group != animation.animation_group) return;

    if (!animation.translation_keyframes.empty()) {
        auto value = vec3f{0, 0, 0};
        switch (animation.interpolation_type) {
            case yocto_interpolation_type::step:
                value = evaluate_keyframed_step(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = evaluate_keyframed_linear(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = evaluate_keyframed_bezier(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            default: log_error("should not have been here");
        }
        for (auto target : animation.node_targets)
            scene.nodes[target].translation = value;
    }
    if (!animation.rotation_keyframes.empty()) {
        auto value = vec4f{0, 0, 0, 1};
        switch (animation.interpolation_type) {
            case yocto_interpolation_type::step:
                value = evaluate_keyframed_step(animation.keyframes_times,
                    animation.rotation_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = evaluate_keyframed_linear(animation.keyframes_times,
                    animation.rotation_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = evaluate_keyframed_bezier(animation.keyframes_times,
                    animation.rotation_keyframes, time);
                break;
        }
        for (auto target : animation.node_targets)
            scene.nodes[target].rotation = value;
    }
    if (!animation.scale_keyframes.empty()) {
        auto value = vec3f{1, 1, 1};
        switch (animation.interpolation_type) {
            case yocto_interpolation_type::step:
                value = evaluate_keyframed_step(
                    animation.keyframes_times, animation.scale_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = evaluate_keyframed_linear(
                    animation.keyframes_times, animation.scale_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = evaluate_keyframed_bezier(
                    animation.keyframes_times, animation.scale_keyframes, time);
                break;
        }
        for (auto target : animation.node_targets)
            scene.nodes[target].scale = value;
    }
}

// Update node transforms
void update_transforms(yocto_scene& scene, yocto_scene_node& node,
    const frame3f& parent = identity_frame3f) {
    auto frame = parent * node.local * make_translation_frame(node.translation) *
                 make_rotation_frame(node.rotation) *
                 make_scaling_frame(node.scale);
    if (node.instance >= 0) scene.instances[node.instance].frame = frame;
    if (node.camera >= 0) scene.cameras[node.camera].frame = frame;
    if (node.environment >= 0)
        scene.environments[node.environment].frame = frame;
    for (auto child : node.children)
        update_transforms(scene, scene.nodes[child], frame);
}

// Update node transforms
void update_transforms(yocto_scene& scene, float time, const string& anim_group) {
    for (auto& agr : scene.animations)
        update_transforms(scene, agr, time, anim_group);
    for (auto& node : scene.nodes) node.children.clear();
    for (auto node_id = 0; node_id < scene.nodes.size(); node_id++) {
        auto& node = scene.nodes[node_id];
        if (node.parent >= 0)
            scene.nodes[node.parent].children.push_back(node_id);
    }
    for (auto& node : scene.nodes)
        if (node.parent >= 0) update_transforms(scene, node);
}

// Compute animation range
vec2f compute_animation_range(const yocto_scene& scene, const string& anim_group) {
    if (scene.animations.empty()) return zero_vec2f;
    auto range = vec2f{+float_max, -float_max};
    for (auto& animation : scene.animations) {
        if (anim_group != "" && animation.animation_group != anim_group)
            continue;
        range[0] = min(range[0], animation.keyframes_times.front());
        range[1] = max(range[1], animation.keyframes_times.back());
    }
    if (range[1] < range[0]) return zero_vec2f;
    return range;
}

// Generate a distribution for sampling a shape uniformly based on area/length.
vector<float> compute_shape_elements_cdf(const yocto_shape& shape) {
    if (!shape.triangles.empty()) {
        return sample_triangles_element_cdf(shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
        return sample_quads_element_cdf(shape.quads, shape.positions);
    } else if (!shape.lines.empty()) {
        return sample_lines_element_cdf(shape.lines, shape.positions);
    } else if (!shape.points.empty()) {
        return sample_points_element_cdf(shape.points.size());
    } else {
        return {};
    }
}

// Sample a shape based on a distribution.
pair<int, vec2f> sample_shape_element(const yocto_shape& shape,
    const vector<float>& elements_cdf, float re, const vec2f& ruv) {
    // TODO: implement sampling without cdf
    if (elements_cdf.empty()) return {};
    if (!shape.triangles.empty()) {
        return sample_triangles_element(elements_cdf, re, ruv);
    } else if (!shape.quads.empty()) {
        return sample_quads_element(elements_cdf, re, ruv);
    } else if (!shape.lines.empty()) {
        return {get<0>(sample_lines_element(elements_cdf, re, ruv[0])), ruv};
    } else if (!shape.points.empty()) {
        return {sample_points_element(elements_cdf, re), ruv};
    } else {
        return {0, zero_vec2f};
    }
}

float sample_shape_element_pdf(const yocto_shape& shape,
    const vector<float>& elements_cdf, int element_id, const vec2f& element_uv) {
    // prob triangle * area triangle = area triangle mesh
    return 1 / elements_cdf.back();
}

// Generate a distribution for sampling a shape uniformly based on area/length.
vector<float> compute_surface_elements_cdf(const yocto_surface& surface) {
    if (!surface.quads_positions.empty()) {
        return sample_quads_element_cdf(
            surface.quads_positions, surface.positions);
    } else {
        return {};
    }
}

// Sample a shape based on a distribution.
pair<int, vec2f> sample_surface_element(const yocto_surface& surface,
    const vector<float>& elements_cdf, float re, const vec2f& ruv) {
    // TODO: implement sampling without cdf
    if (elements_cdf.empty()) return {};
    if (!surface.quads_positions.empty()) {
        return sample_quads_element(elements_cdf, re, ruv);
    } else {
        return {0, zero_vec2f};
    }
}

float sample_surface_element_pdf(const yocto_surface& surface,
    const vector<float>& elements_cdf, int element_id, const vec2f& element_uv) {
    // prob triangle * area triangle = area triangle mesh
    return 1 / elements_cdf.back();
}

// Update environment CDF for sampling.
vector<float> compute_environment_texels_cdf(
    const yocto_scene& scene, const yocto_environment& environment) {
    if (environment.emission_texture < 0) return {};
    auto& texture  = scene.textures[environment.emission_texture];
    auto  size     = evaluate_texture_size(texture);
    auto  elem_cdf = vector<float>(size[0] * size[1]);
    if (size != zero_vec2i) {
        for (auto i = 0; i < elem_cdf.size(); i++) {
            auto ij     = vec2i{i % size[0], i / size[0]};
            auto th     = (ij[1] + 0.5f) * pif / size[1];
            auto value  = lookup_texture(texture, ij);
            elem_cdf[i] = max(make_shorter_vec(value)) * sin(th);
            if (i) elem_cdf[i] += elem_cdf[i - 1];
        }
    } else {
        log_error("empty texture");
    }
    return elem_cdf;
}

// Sample an environment based on texels
vec3f sample_environment_direction(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    float re, const vec2f& ruv) {
    if (!texels_cdf.empty() && environment.emission_texture >= 0) {
        auto& texture = scene.textures[environment.emission_texture];
        auto  idx     = sample_discrete_distribution(texels_cdf, re);
        auto  size    = evaluate_texture_size(texture);
        auto  u       = (idx % size[0] + 0.5f) / size[0];
        auto  v       = (idx / size[0] + 0.5f) / size[1];
        return evaluate_environment_direction(environment, {u, v});
    } else {
        return sample_sphere_direction(ruv);
    }
}

// Sample an environment based on texels
float sample_environment_direction_pdf(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    const vec3f& direction) {
    if (!texels_cdf.empty() && environment.emission_texture >= 0) {
        auto& texture = scene.textures[environment.emission_texture];
        auto  size    = evaluate_texture_size(texture);
        auto texcoord = evaluate_environment_texturecoord(environment, direction);
        auto i        = (int)(texcoord[0] * size[0]);
        auto j        = (int)(texcoord[1] * size[1]);
        auto idx      = j * size[0] + i;
        auto prob     = sample_discrete_distribution_pdf(texels_cdf, idx) /
                    texels_cdf.back();
        auto angle = (2 * pif / size[0]) * (pif / size[1]) *
                     sin(pif * (j + 0.5f) / size[1]);
        return prob / angle;
    } else {
        return sample_sphere_direction_pdf(direction);
    }
}

// Build a scene BVH
bvh_scene make_scene_bvh(
    const yocto_scene& scene, const build_bvh_options& options) {
    auto scope = log_trace_scoped("building scene bvh");
    // shapes
    auto shape_bvhs = vector<bvh_shape>();
    for (auto& shape : scene.shapes) {
        auto shape_bvh = bvh_shape{};
        if (!shape.points.empty()) {
            shape_bvh = make_shape_bvh(
                shape.points, shape.positions, shape.radius);
        } else if (!shape.lines.empty()) {
            shape_bvh = make_shape_bvh(
                shape.lines, shape.positions, shape.radius);
        } else if (!shape.triangles.empty()) {
            shape_bvh = make_shape_bvh(shape.triangles, shape.positions);
        } else if (!shape.quads.empty()) {
            shape_bvh = make_shape_bvh(shape.quads, shape.positions);
        } else {
            shape_bvh = {};
        }
        shape_bvhs.push_back(shape_bvh);
    }

    // surfaces
    auto surface_bvhs = vector<bvh_shape>();
    for (auto& surface : scene.surfaces) {
        auto surface_bvh = bvh_shape{};
        if (!surface.quads_positions.empty()) {
            surface_bvh = make_shape_bvh(
                surface.quads_positions, surface.positions);
        } else {
            surface_bvh = {};
        }
        surface_bvhs.push_back(surface_bvh);
    }

    // instances
    auto bvh_instances = vector<bvh_instance>{};
    for (auto& instance : scene.instances) {
        bvh_instances.push_back({instance.frame, inverse(instance.frame, false),
            instance.shape, instance.surface});
    }

    // build bvh
    auto bvh = make_scene_bvh(bvh_instances, shape_bvhs, surface_bvhs);
    build_scene_bvh(bvh, options);
    return bvh;
}

// Refits a scene BVH
void refit_scene_bvh(const yocto_scene& scene, bvh_scene& bvh,
    const vector<int>& updated_instances, const vector<int>& updated_shapes,
    const vector<int>& updated_surfaces) {
    for (auto shape_id : updated_shapes)
        update_shape_bvh(get_shape_bvh(bvh, shape_id),
            scene.shapes[shape_id].positions, scene.shapes[shape_id].radius);
    for (auto surface_id : updated_surfaces)
        update_shape_bvh(get_surface_bvh(bvh, surface_id),
            scene.surfaces[surface_id].positions);

    auto bvh_instances = vector<bvh_instance>{};
    for (auto& instance : scene.instances) {
        bvh_instances.push_back({instance.frame, inverse(instance.frame),
            instance.shape, instance.surface});
    }
    update_scene_bvh(bvh, bvh_instances);

    refit_scene_bvh(bvh, updated_instances, updated_shapes, updated_surfaces);
}

// Add missing names and resolve duplicated names.
void add_missing_names(yocto_scene& scene) {
    auto fix_names = [](auto& vals, const string& base) {
        auto nmap = unordered_map<string, int>();
        for (auto& value : vals) {
            if (value.name == "") value.name = base;
            if (nmap.find(value.name) == nmap.end()) {
                nmap[value.name] = 0;
            } else {
                nmap[value.name] += 1;
                value.name = value.name + "_" + std::to_string(nmap[value.name]);
            }
        }
    };
    fix_names(scene.cameras, "camera");
    fix_names(scene.shapes, "shape");
    fix_names(scene.surfaces, "surface");
    fix_names(scene.textures, "texture");
    fix_names(scene.voltextures, "voltexture");
    fix_names(scene.materials, "material");
    fix_names(scene.instances, "instance");
    fix_names(scene.environments, "environment");
    fix_names(scene.nodes, "node");
    fix_names(scene.animations, "animation");
}

// Add missing tangent space if needed.
void add_missing_tangent_space(yocto_scene& scene) {
    for (auto& shape : scene.shapes) {
        auto& material = scene.materials[shape.material];
        if (!shape.tangentspaces.empty() || shape.texturecoords.empty())
            continue;
        if (material.normal_texture < 0 && material.bump_texture < 0) continue;
        if (!shape.triangles.empty()) {
            if (shape.normals.empty())
                shape.normals = compute_vertex_normals(
                    shape.triangles, shape.positions);
            shape.tangentspaces = compute_tangent_spaces(shape.triangles,
                shape.positions, shape.normals, shape.texturecoords);
        } else {
            log_error("type not supported");
        }
    }
}

// Add missing materials.
void add_missing_materials(yocto_scene& scene) {
    auto material_id = -1;
    for (auto& shape : scene.shapes) {
        if (shape.material >= 0) continue;
        if (material_id < 0) {
            auto material    = yocto_material{};
            material.name    = "<default>";
            material.diffuse = {0.2f, 0.2f, 0.2f};
            scene.materials.push_back(material);
            material_id = (int)scene.materials.size() - 1;
        }
        shape.material = material_id;
    }
    for (auto& surface : scene.surfaces) {
        if (!surface.materials.empty()) continue;
        if (material_id < 0) {
            auto material    = yocto_material{};
            material.name    = "<default>";
            material.diffuse = {0.2f, 0.2f, 0.2f};
            scene.materials.push_back(material);
            material_id = (int)scene.materials.size() - 1;
        }
        surface.materials.push_back(material_id);
    }
}

// Add missing cameras.
void add_missing_cameras(yocto_scene& scene) {
    if (scene.cameras.empty()) {
        auto camera = yocto_camera{};
        camera.name = "<view>";
        set_camera_view(camera, compute_scene_bounds(scene), {0, 0, 1});
        scene.cameras.push_back(camera);
    }
}

// Add a sky environment
void add_sky_environment(yocto_scene& scene, float sun_angle) {
    auto texture      = yocto_texture{};
    texture.name      = "<sky>";
    texture.filename  = "textures/sky.hdr";
    texture.hdr_image = make_sunsky_image({1024, 512}, sun_angle);
    scene.textures.push_back(texture);
    auto environment             = yocto_environment{};
    environment.name             = "<sky>";
    environment.emission         = {1, 1, 1};
    environment.emission_texture = (int)scene.textures.size() - 1;
    scene.environments.push_back(environment);
}

// Checks for validity of the scene.
vector<string> validate_scene(const yocto_scene& scene, bool skip_textures) {
    auto errs        = vector<string>();
    auto check_names = [&errs](const auto& vals, const string& base) {
        auto used = unordered_map<string, int>();
        for (auto& value : vals) used[value.name] += 1;
        for (auto& [name, used] : used) {
            if (name == "") {
                errs.push_back("empty " + base + " name");
            } else if (used > 1) {
                errs.push_back("duplicated " + base + " name " + name);
            }
        }
    };
    auto check_empty_textures = [&errs](const vector<yocto_texture>& vals) {
        for (auto& value : vals) {
            if (value.hdr_image.empty() && value.ldr_image.empty()) {
                errs.push_back("empty texture " + value.name);
            }
        }
    };

    check_names(scene.cameras, "camera");
    check_names(scene.shapes, "shape");
    check_names(scene.surfaces, "surface");
    check_names(scene.textures, "texture");
    check_names(scene.voltextures, "voltexture");
    check_names(scene.materials, "material");
    check_names(scene.instances, "instance");
    check_names(scene.environments, "environment");
    check_names(scene.nodes, "node");
    check_names(scene.animations, "animation");
    if (!skip_textures) check_empty_textures(scene.textures);

    return errs;
}

// Logs validations errors
void log_validation_errors(const yocto_scene& scene, bool skip_textures) {
    for (auto err : validate_scene(scene, skip_textures))
        log_error(err + " [validation]");
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Scene intersection.
scene_intersection intersect_scene(const yocto_scene& scene,
    const bvh_scene& bvh, const ray3f& ray, bool find_any) {
    auto isec = scene_intersection();
    if (!intersect_scene_bvh(bvh, ray, find_any, isec.distance,
            isec.instance_id, isec.element_id, isec.element_uv))
        return {};
    return isec;
}

// Instance intersection.
scene_intersection intersect_scene(const yocto_scene& scene, int instance_id,
    const bvh_scene& bvh, const ray3f& ray, bool find_any) {
    auto& instance = scene.instances[instance_id];
    auto  isec     = scene_intersection();
    auto  tray     = transform_ray_inverse(instance.frame, ray);
    if (instance.shape >= 0) {
        if (!intersect_shape_bvh(bvh.shape_bvhs[instance.shape], tray, find_any,
                isec.distance, isec.element_id, isec.element_uv))
            return {};
    } else if (instance.surface) {
        if (!intersect_shape_bvh(bvh.shape_bvhs[instance.surface], tray,
                find_any, isec.distance, isec.element_id, isec.element_uv))
            return {};
    } else {
        return {};
    }
    isec.instance_id = instance_id;
    return isec;
}

// Shape element normal.
vec3f evaluate_shape_element_normal(const yocto_shape& shape, int element_id) {
    auto norm = zero_vec3f;
    if (!shape.triangles.empty()) {
        auto t = shape.triangles[element_id];
        norm   = triangle_normal(shape.positions[t[0]], shape.positions[t[1]],
            shape.positions[t[2]]);
    } else if (!shape.quads.empty()) {
        auto q = shape.quads[element_id];
        norm   = quad_normal(shape.positions[q[0]], shape.positions[q[1]],
            shape.positions[q[2]], shape.positions[q[3]]);
    } else if (!shape.lines.empty()) {
        auto l = shape.lines[element_id];
        norm   = line_tangent(shape.positions[l[0]], shape.positions[l[1]]);
    } else {
        norm = {0, 0, 1};
    }
    return norm;
}

// Shape element normal.
vec4f evaluate_shape_element_tangentspace(
    const yocto_shape& shape, int element_id) {
    auto tangsp = zero_vec4f;
    if (!shape.triangles.empty()) {
        auto t    = shape.triangles[element_id];
        auto norm = triangle_normal(shape.positions[t[0]],
            shape.positions[t[1]], shape.positions[t[2]]);
        auto txty = pair<vec3f, vec3f>();
        if (shape.texturecoords.empty()) {
            txty = triangle_tangents_fromuv(shape.positions[t[0]],
                shape.positions[t[1]], shape.positions[t[2]], {0, 0}, {1, 0},
                {0, 1});
        } else {
            txty = triangle_tangents_fromuv(shape.positions[t[0]],
                shape.positions[t[1]], shape.positions[t[2]],
                shape.texturecoords[t[0]], shape.texturecoords[t[1]],
                shape.texturecoords[t[2]]);
        }
        auto tx = txty.first, ty = txty.second;
        tx     = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        tangsp = {tx[0], tx[1], tx[2], s};
    }
    return tangsp;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T evaluate_shape_elem(const yocto_shape& shape, const vector<T>& vals,
    int element_id, const vec2f& element_uv) {
    if (vals.empty()) return {};
    if (!shape.triangles.empty()) {
        auto t = shape.triangles[element_id];
        return interpolate_triangle(
            vals[t[0]], vals[t[1]], vals[t[2]], element_uv);
    } else if (!shape.quads.empty()) {
        auto q = shape.quads[element_id];
        if (q[3] == q[2])
            return interpolate_triangle(
                vals[q[0]], vals[q[1]], vals[q[2]], element_uv);
        return interpolate_quad(
            vals[q[0]], vals[q[1]], vals[q[2]], vals[q[3]], element_uv);
    } else if (!shape.lines.empty()) {
        auto l = shape.lines[element_id];
        return interpolate_line(vals[l[0]], vals[l[1]], element_uv[0]);
    } else if (!shape.points.empty()) {
        return vals[shape.points[element_id]];
    } else {
        return {};
    }
}

// Shape values interpolated using barycentric coordinates
vec3f evaluate_shape_position(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    return evaluate_shape_elem(shape, shape.positions, element_id, element_uv);
}
vec3f evaluate_shape_normal(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.normals.empty())
        return evaluate_shape_element_normal(shape, element_id);
    return normalize(
        evaluate_shape_elem(shape, shape.normals, element_id, element_uv));
}
vec2f evaluate_shape_texturecoord(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.texturecoords.empty()) return element_uv;
    return evaluate_shape_elem(
        shape, shape.texturecoords, element_id, element_uv);
}
vec4f evaluate_shape_color(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.colors.empty()) return {1, 1, 1, 1};
    return evaluate_shape_elem(shape, shape.colors, element_id, element_uv);
}
float evaluate_shape_radius(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.radius.empty()) return 0.001f;
    return evaluate_shape_elem(shape, shape.radius, element_id, element_uv);
}
vec4f evaluate_shape_tangentspace(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.tangentspaces.empty())
        return evaluate_shape_element_tangentspace(shape, element_id);
    return evaluate_shape_elem(
        shape, shape.tangentspaces, element_id, element_uv);
}
vec3f evaluate_shape_tangentspace(const yocto_shape& shape, int element_id,
    const vec2f& element_uv, bool& left_handed) {
    auto tangsp = (shape.tangentspaces.empty()) ?
                      evaluate_shape_element_tangentspace(shape, element_id) :
                      evaluate_shape_elem(
                          shape, shape.tangentspaces, element_id, element_uv);
    left_handed = tangsp[3] < 0;
    return {tangsp[0], tangsp[1], tangsp[2]};
}
// Shading normals including material perturbations.
vec3f evaluate_shape_shading_normal(const yocto_scene& scene,
    const yocto_shape& shape, int element_id, const vec2f& element_uv,
    const vec3f& outgoing) {
    if (!shape.triangles.empty()) {
        auto  normal   = evaluate_shape_normal(shape, element_id, element_uv);
        auto& material = scene.materials[shape.material];
        if (material.normal_texture >= 0) {
            auto texcoord = evaluate_shape_texturecoord(
                shape, element_id, element_uv);
            auto  left_handed    = false;
            auto& normal_texture = scene.textures[material.normal_texture];
            auto  texture        = make_shorter_vec(
                evaluate_texture(normal_texture, texcoord));
            texture    = texture * 2 - vec3f{1, 1, 1};
            texture[1] = -texture[1];  // flip vertical axis to align green with
                                       // image up
            auto tu = orthonormalize(evaluate_shape_tangentspace(shape,
                                         element_id, element_uv, left_handed),
                normal);
            auto tv = normalize(cross(normal, tu) * (left_handed ? -1.0f : 1.0f));
            normal = normalize(
                texture[0] * tu + texture[1] * tv + texture[2] * normal);
        }
        if (material.double_sided && dot(normal, outgoing) < 0)
            normal = -normal;
        return normal;
    } else if (!shape.quads.empty()) {
        auto  normal   = evaluate_shape_normal(shape, element_id, element_uv);
        auto& material = scene.materials[shape.material];
        if (material.double_sided && dot(normal, outgoing) < 0)
            normal = -normal;
        return normal;
    } else if (!shape.lines.empty()) {
        return orthonormalize(
            outgoing, evaluate_shape_normal(shape, element_id, element_uv));
    } else {
        return outgoing;
    }
}

// Shape element normal.
vec3f evaluate_surface_element_normal(
    const yocto_surface& surface, int element_id) {
    auto norm = zero_vec3f;
    if (!surface.quads_positions.empty()) {
        auto q = surface.quads_positions[element_id];
        norm   = quad_normal(surface.positions[q[0]], surface.positions[q[1]],
            surface.positions[q[2]], surface.positions[q[3]]);
    } else {
        norm = {0, 0, 1};
    }
    return norm;
}

// Shape element normal.
vec4f evaluate_surface_element_tangentspace(
    const yocto_surface& surface, int element_id) {
    return zero_vec4f;
}

// override for face-varying data
template <typename T>
T evaluate_surface_elem(const yocto_surface& surface, const vector<T>& vals,
    const vector<vec4i>& quads, int element_id, const vec2f& element_uv) {
    if (vals.empty()) return {};
    auto q = quads[element_id];
    if (q[3] == q[2])
        return interpolate_triangle(
            vals[q[0]], vals[q[1]], vals[q[2]], element_uv);
    return interpolate_quad(
        vals[q[0]], vals[q[1]], vals[q[2]], vals[q[3]], element_uv);
}

// Shape values interpolated using barycentric coordinates
vec3f evaluate_surface_position(
    const yocto_surface& surface, int element_id, const vec2f& element_uv) {
    return evaluate_surface_elem(surface, surface.positions,
        surface.quads_positions, element_id, element_uv);
}
vec3f evaluate_surface_normal(
    const yocto_surface& surface, int element_id, const vec2f& element_uv) {
    if (surface.normals.empty())
        return evaluate_surface_element_normal(surface, element_id);
    return normalize(evaluate_surface_elem(surface, surface.normals,
        surface.quads_normals, element_id, element_uv));
}
vec2f evaluate_surface_texturecoord(
    const yocto_surface& surface, int element_id, const vec2f& element_uv) {
    if (surface.texturecoords.empty()) return element_uv;
    return evaluate_surface_elem(surface, surface.texturecoords,
        surface.quads_texturecoords, element_id, element_uv);
}
vec4f evaluate_surface_tangentspace(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    return evaluate_shape_element_tangentspace(shape, element_id);
}
vec3f evaluate_surface_tangentspace(const yocto_surface& surface,
    int element_id, const vec2f& element_uv, bool& left_handed) {
    auto tangsp = evaluate_surface_element_tangentspace(surface, element_id);
    left_handed = tangsp[3] < 0;
    return {tangsp[0], tangsp[1], tangsp[2]};
}
// Shading normals including material perturbations.
vec3f evaluate_surface_shading_normal(const yocto_scene& scene,
    const yocto_surface& surface, int element_id, const vec2f& element_uv,
    const vec3f& outgoing) {
    if (!surface.quads_positions.empty()) {
        auto  normal = evaluate_surface_normal(surface, element_id, element_uv);
        auto& material = scene.materials[get_surface_element_material(
            surface, element_id)];
        if (material.double_sided && dot(normal, outgoing) < 0)
            normal = -normal;
        return normal;
    } else {
        return outgoing;
    }
}
// Per-element material.
int get_surface_element_material(const yocto_surface& surface, int element_id) {
    if (surface.materials.empty()) return -1;
    if (surface.quads_materials.empty()) return surface.materials.front();
    return surface.materials[surface.quads_materials[element_id]];
}

// Instance values interpolated using barycentric coordinates.
vec3f evaluate_instance_position(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        return transform_point(instance.frame,
            evaluate_shape_position(
                scene.shapes[instance.shape], element_id, element_uv));
    } else if (instance.surface >= 0) {
        return transform_point(instance.frame,
            evaluate_surface_position(
                scene.surfaces[instance.surface], element_id, element_uv));
    } else {
        log_error("empty instance");
        return zero_vec3f;
    }
}
vec3f evaluate_instance_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        return transform_direction(
            instance.frame, evaluate_shape_normal(scene.shapes[instance.shape],
                                element_id, element_uv));
    } else if (instance.surface >= 0) {
        return transform_direction(instance.frame,
            evaluate_surface_normal(
                scene.surfaces[instance.surface], element_id, element_uv));
    } else {
        log_error("empty instance");
        return zero_vec3f;
    }
}
vec2f evaluate_instance_texturecoord(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        return evaluate_shape_texturecoord(
            scene.shapes[instance.shape], element_id, element_uv);
    } else if (instance.surface >= 0) {
        return evaluate_surface_texturecoord(
            scene.surfaces[instance.surface], element_id, element_uv);
    } else {
        log_error("empty instance");
        return zero_vec2f;
    }
}
vec3f evaluate_instance_tangentspace(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    bool& left_handed) {
    if (instance.shape >= 0) {
        return transform_direction(instance.frame,
            evaluate_shape_tangentspace(scene.shapes[instance.shape],
                element_id, element_uv, left_handed));
    } else if (instance.surface >= 0) {
        return transform_direction(instance.frame,
            evaluate_surface_tangentspace(scene.surfaces[instance.surface],
                element_id, element_uv, left_handed));
    } else {
        log_error("empty instance");
        return zero_vec3f;
    }
}
// Instance element values.
vec3f evaluate_instance_element_normal(
    const yocto_scene& scene, const yocto_instance& instance, int element_id) {
    if (instance.shape >= 0) {
        return transform_direction(
            instance.frame, evaluate_shape_element_normal(
                                scene.shapes[instance.shape], element_id));
    } else if (instance.surface >= 0) {
        return transform_direction(
            instance.frame, evaluate_surface_element_normal(
                                scene.surfaces[instance.surface], element_id));
    } else {
        log_error("empty instance");
        return zero_vec3f;
    }
}
// Shading normals including material perturbations.
vec3f evaluate_instance_shading_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    const vec3f& outgoing) {
    if (instance.shape >= 0) {
        return transform_direction(instance.frame,
            evaluate_shape_shading_normal(scene, scene.shapes[instance.shape],
                element_id, element_uv,
                transform_direction_inverse(instance.frame, outgoing)));
    } else if (instance.surface >= 0) {
        return transform_direction(instance.frame,
            evaluate_surface_shading_normal(scene,
                scene.surfaces[instance.surface], element_id, element_uv,
                transform_direction_inverse(instance.frame, outgoing)));
    } else {
        log_error("empty instance");
        return zero_vec3f;
    }
}

// Material values
vec3f evaluate_instance_emission(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        auto& shape = scene.shapes[instance.shape];
        return evaluate_material_emission(scene, scene.materials[shape.material],
            evaluate_shape_texturecoord(shape, element_id, element_uv),
            evaluate_shape_color(shape, element_id, element_uv));
    } else if (instance.surface >= 0) {
        auto& surface = scene.surfaces[instance.surface];
        return evaluate_material_emission(scene,
            scene.materials[get_surface_element_material(surface, element_id)],
            evaluate_surface_texturecoord(surface, element_id, element_uv),
            {1, 1, 1, 1});
    } else {
        return zero_vec3f;
    }
}
microfacet_brdf evaluate_instance_brdf(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        auto& shape = scene.shapes[instance.shape];
        return evaluate_material_brdf(scene, scene.materials[shape.material],
            evaluate_shape_texturecoord(shape, element_id, element_uv),
            evaluate_shape_color(shape, element_id, element_uv));
    } else if (instance.surface >= 0) {
        auto& surface = scene.surfaces[instance.surface];
        return evaluate_material_brdf(scene,
            scene.materials[get_surface_element_material(surface, element_id)],
            evaluate_surface_texturecoord(surface, element_id, element_uv),
            {1, 1, 1, 1});
    } else {
        return {};
    }
}
float evaluate_instance_opacity(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    if (instance.shape >= 0) {
        auto& shape = scene.shapes[instance.shape];
        return evaluate_material_opacity(scene, scene.materials[shape.material],
            evaluate_shape_texturecoord(shape, element_id, element_uv),
            evaluate_shape_color(shape, element_id, element_uv));
    } else if (instance.surface >= 0) {
        auto& surface = scene.surfaces[instance.surface];
        return evaluate_material_opacity(scene,
            scene.materials[get_surface_element_material(surface, element_id)],
            evaluate_surface_texturecoord(surface, element_id, element_uv),
            {1, 1, 1, 1});
    } else {
        return 0;
    }
}
bool is_instance_emissive(
    const yocto_scene& scene, const yocto_instance& instance) {
    if (instance.shape >= 0) {
        auto& shape = scene.shapes[instance.shape];
        return scene.materials[shape.material].emission != zero_vec3f;
    } else if (instance.surface >= 0) {
        auto& surface = scene.surfaces[instance.surface];
        for (auto material_id : surface.materials)
            if (scene.materials[material_id].emission != zero_vec3f)
                return true;
        return false;
    } else {
        return false;
    }
}

// Environment texture coordinates from the direction.
vec2f evaluate_environment_texturecoord(
    const yocto_environment& environment, const vec3f& direction) {
    auto wl = transform_direction_inverse(environment.frame, direction);
    auto environment_uv = vec2f{
        atan2(wl[2], wl[0]) / (2 * pif), acos(clamp(wl[1], -1.0f, 1.0f)) / pif};
    if (environment_uv[0] < 0) environment_uv[0] += 1;
    return environment_uv;
}
// Evaluate the environment direction.
vec3f evaluate_environment_direction(
    const yocto_environment& environment, const vec2f& environment_uv) {
    return transform_direction(environment.frame,
        {cos(environment_uv[0] * 2 * pif) * sin(environment_uv[1] * pif),
            cos(environment_uv[1] * pif),
            sin(environment_uv[0] * 2 * pif) * sin(environment_uv[1] * pif)});
}
// Evaluate the environment color.
vec3f evaluate_environment_emission(const yocto_scene& scene,
    const yocto_environment& environment, const vec3f& direction) {
    auto ke = environment.emission;
    if (environment.emission_texture >= 0) {
        auto& emission_texture = scene.textures[environment.emission_texture];
        ke *= make_shorter_vec(evaluate_texture(emission_texture,
            evaluate_environment_texturecoord(environment, direction)));
    }
    return ke;
}
// Evaluate all environment color.
vec3f evaluate_environment_emission(
    const yocto_scene& scene, const vec3f& direction) {
    auto ke = zero_vec3f;
    for (auto& environment : scene.environments)
        ke += evaluate_environment_emission(scene, environment, direction);
    return ke;
}

// Check texture size
vec2i evaluate_texture_size(const yocto_texture& texture) {
    if (!texture.hdr_image.empty()) {
        return texture.hdr_image.size();
    } else if (!texture.ldr_image.empty()) {
        return texture.ldr_image.size();
    } else {
        return zero_vec2i;
    }
}

// Lookup a texture value
vec4f lookup_texture(const yocto_texture& texture, const vec2i& ij) {
    if (!texture.hdr_image.empty()) {
        return texture.hdr_image[ij];
    } else if (!texture.ldr_image.empty() && !texture.ldr_as_linear) {
        return srgb_to_linear(byte_to_float(texture.ldr_image[ij]));
    } else if (!texture.ldr_image.empty() && texture.ldr_as_linear) {
        return byte_to_float(texture.ldr_image[ij]);
    } else {
        return zero_vec4f;
    }
}

// Evaluate a texture
vec4f evaluate_texture(const yocto_texture& texture, const vec2f& texcoord) {
    if (texture.hdr_image.empty() && texture.ldr_image.empty())
        return {1, 1, 1, 1};

    // get image width/height
    auto size  = evaluate_texture_size(texture);
    auto width = size[0], height = size[1];

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (texture.clamp_to_edge) {
        s = clamp(texcoord[0], 0.0f, 1.0f) * width;
        t = clamp(texcoord[1], 0.0f, 1.0f) * height;
    } else {
        s = fmod(texcoord[0], 1.0f) * width;
        if (s < 0) s += width;
        t = fmod(texcoord[1], 1.0f) * height;
        if (t < 0) t += height;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, width - 1), j = clamp((int)t, 0, height - 1);
    auto ii = (i + 1) % width, jj = (j + 1) % height;
    auto u = s - i, v = t - j;

    // nearest-neighbor interpolation
    if (texture.no_interpolation) {
        i = u < 0.5 ? i : min(i + 1, width - 1);
        j = v < 0.5 ? j : min(j + 1, height - 1);
        return lookup_texture(texture, {i, j});
    }

    // handle interpolation
    return lookup_texture(texture, {i, j}) * (1 - u) * (1 - v) +
           lookup_texture(texture, {i, jj}) * (1 - u) * v +
           lookup_texture(texture, {ii, j}) * u * (1 - v) +
           lookup_texture(texture, {ii, jj}) * u * v;
}

// Lookup a texture value
float lookup_voltexture(const yocto_voltexture& texture, const vec3i& ijk) {
    if (!texture.volume_data.empty()) {
        return texture.volume_data[ijk];
    } else {
        return 0;
    }
}

// Evaluate a volume texture
float evaluate_voltexture(const yocto_voltexture& texture, const vec3f& texcoord) {
    if (texture.volume_data.empty()) return 1;

    // get image width/height
    auto width  = texture.volume_data.width();
    auto height = texture.volume_data.height();
    auto depth  = texture.volume_data.depth();

    // get coordinates normalized for tiling
    auto s = clamp((texcoord[0] + 1.0f) * 0.5f, 0.0f, 1.0f) * width;
    auto t = clamp((texcoord[1] + 1.0f) * 0.5f, 0.0f, 1.0f) * height;
    auto r = clamp((texcoord[2] + 1.0f) * 0.5f, 0.0f, 1.0f) * depth;

    // get image coordinates and residuals
    auto i  = clamp((int)s, 0, width - 1);
    auto j  = clamp((int)t, 0, height - 1);
    auto k  = clamp((int)r, 0, depth - 1);
    auto ii = (i + 1) % width, jj = (j + 1) % height, kk = (k + 1) % depth;
    auto u = s - i, v = t - j, w = r - k;

    // nearest-neighbor interpolation
    if (texture.no_interpolation) {
        i = u < 0.5 ? i : min(i + 1, width - 1);
        j = v < 0.5 ? j : min(j + 1, height - 1);
        k = w < 0.5 ? k : min(k + 1, depth - 1);
        return lookup_voltexture(texture, {i, j, k});
    }

    // trilinear interpolation
    return lookup_voltexture(texture, {i, j, k}) * (1 - u) * (1 - v) * (1 - w) +
           lookup_voltexture(texture, {ii, j, k}) * u * (1 - v) * (1 - w) +
           lookup_voltexture(texture, {i, jj, k}) * (1 - u) * v * (1 - w) +
           lookup_voltexture(texture, {i, j, kk}) * (1 - u) * (1 - v) * w +
           lookup_voltexture(texture, {i, jj, kk}) * (1 - u) * v * w +
           lookup_voltexture(texture, {ii, j, kk}) * u * (1 - v) * w +
           lookup_voltexture(texture, {ii, jj, k}) * u * v * (1 - w) +
           lookup_voltexture(texture, {ii, jj, kk}) * u * v * w;
}

// Set and evaluate camera parameters. Setters take zeros as default values.
float get_camera_fovx(const yocto_camera& camera) {
    return 2 * atan(camera.film_size[0] / (2 * camera.focal_length));
}
float get_camera_fovy(const yocto_camera& camera) {
    return 2 * atan(camera.film_size[1] / (2 * camera.focal_length));
}
float get_camera_aspect(const yocto_camera& camera) {
    return camera.film_size[0] / camera.film_size[1];
}
vec2i get_camera_image_size(const yocto_camera& camera, const vec2i& size) {
    return get_image_size(size, camera.film_size[0] / camera.film_size[1]);
}
void set_camera_fovy(yocto_camera& camera, float fovy, float aspect, float width) {
    camera.film_size    = {width, width / aspect};
    camera.focal_length = camera.film_size[1] / (2 * tan(fovy / 2));
}

// add missing camera
void set_camera_view(yocto_camera& camera, const bbox3f& bbox,
    const vec3f& view_direction, const vec2f& film, float focal) {
    camera.orthographic = false;
    if (film != zero_vec2f) camera.film_size = film;
    if (focal != 0) camera.focal_length = focal;
    auto bbox_center = (bbox.max + bbox.min) / 2.0f;
    auto bbox_radius = length(bbox.max - bbox.min) / 2;
    auto camera_dir  = (view_direction == zero_vec3f) ?
                          camera.frame.origin - bbox_center :
                          view_direction;
    if (camera_dir == zero_vec3f) camera_dir = {0, 0, 1};
    auto camera_fov = min(get_camera_fovx(camera), get_camera_fovy(camera));
    if (camera_fov == 0) camera_fov = 45 * pif / 180;
    auto camera_dist      = bbox_radius / sin(camera_fov / 2);
    auto from             = camera_dir * (camera_dist * 1) + bbox_center;
    auto to               = bbox_center;
    auto up               = vec3f{0, 1, 0};
    camera.frame          = make_lookat_frame(from, to, up);
    camera.focus_distance = length(from - to);
    camera.lens_aperture  = 0;
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f evaluate_camera_ray(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
    auto distance = camera.focal_length;
    if (camera.focus_distance < float_max) {
        distance = camera.focal_length * camera.focus_distance /
                   (camera.focus_distance - camera.focal_length);
    }
    auto e   = vec3f{lens_uv[0] * camera.lens_aperture,
        lens_uv[1] * camera.lens_aperture, 0};
    auto q   = vec3f{camera.film_size[0] * (0.5f - image_uv[0]),
        camera.film_size[1] * (image_uv[1] - 0.5f), distance};
    auto ray = make_ray(transform_point(camera.frame, e),
        transform_direction(camera.frame, normalize(e - q)));
    return ray;
}

// Generates a ray from a camera.
ray3f evaluate_camera_ray(const yocto_camera& camera, const vec2i& image_ij,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv) {
    auto image_uv = vec2f{(image_ij[0] + pixel_uv[0]) / image_size[0],
        (image_ij[1] + pixel_uv[1]) / image_size[1]};
    return evaluate_camera_ray(camera, image_uv, lens_uv);
}

// Generates a ray from a camera.
ray3f evaluate_camera_ray(const yocto_camera& camera, int idx,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv) {
    auto image_ij = vec2i{idx % image_size[0], idx / image_size[0]};
    auto image_uv = vec2f{(image_ij[0] + pixel_uv[0]) / image_size[0],
        (image_ij[1] + pixel_uv[1]) / image_size[1]};
    return evaluate_camera_ray(camera, image_uv, lens_uv);
}

// Evaluates material parameters.
vec3f evaluate_material_emission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    auto emission = material.emission * make_shorter_vec(shape_color);
    if (material.emission_texture >= 0) {
        auto& emission_texture = scene.textures[material.emission_texture];
        emission *= make_shorter_vec(
            evaluate_texture(emission_texture, texturecoord));
    }
    return emission;
}
vec3f evaluate_material_diffuse(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    if (!material.base_metallic) {
        auto diffuse = material.diffuse * make_shorter_vec(shape_color);
        if (material.diffuse_texture >= 0) {
            auto& diffuse_texture = scene.textures[material.diffuse_texture];
            diffuse *= make_shorter_vec(
                evaluate_texture(diffuse_texture, texturecoord));
        }
        return diffuse;
    } else {
        auto base = material.diffuse * make_shorter_vec(shape_color);
        if (material.diffuse_texture >= 0) {
            auto& diffuse_texture = scene.textures[material.diffuse_texture];
            base *= make_shorter_vec(
                evaluate_texture(diffuse_texture, texturecoord));
        }
        auto metallic = material.specular;
        if (material.specular_texture >= 0) {
            auto& specular_texture = scene.textures[material.specular_texture];
            metallic *= evaluate_texture(specular_texture, texturecoord)[2];
        }
        return base * (1 - metallic);
    }
}
vec3f evaluate_material_specular(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    if (!material.base_metallic) {
        auto specular = material.specular * make_shorter_vec(shape_color);
        if (material.specular_texture >= 0) {
            auto& specular_texture = scene.textures[material.specular_texture];
            specular *= make_shorter_vec(
                evaluate_texture(specular_texture, texturecoord));
        }
        return specular;
    } else {
        auto base = material.diffuse * make_shorter_vec(shape_color);
        if (material.diffuse_texture >= 0) {
            auto& diffuse_texture = scene.textures[material.diffuse_texture];
            base *= make_shorter_vec(
                evaluate_texture(diffuse_texture, texturecoord));
        }
        auto metallic = material.specular[0];
        if (material.specular_texture >= 0) {
            auto& specular_texture = scene.textures[material.specular_texture];
            metallic *= evaluate_texture(specular_texture, texturecoord)[2];
        }
        return base * metallic + 0.04f * (1 - metallic);
    }
}
float evaluate_material_roughness(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    if (!material.base_metallic) {
        if (!material.gltf_textures) {
            auto roughness = material.roughness;
            if (material.roughness_texture >= 0) {
                auto& roughness_texture = scene.textures[material.roughness_texture];
                roughness *= evaluate_texture(roughness_texture, texturecoord)[0];
            }
            return roughness * roughness;
        } else {
            auto glossiness = 1 - material.roughness;
            if (material.roughness_texture >= 0) {
                auto& roughness_texture = scene.textures[material.roughness_texture];
                glossiness *= evaluate_texture(
                    roughness_texture, texturecoord)[3];
            }
            auto roughness = 1 - glossiness;
            return roughness * roughness;
        }
    } else {
        auto roughness = material.roughness;
        if (material.roughness_texture >= 0) {
            auto& roughness_texture = scene.textures[material.roughness_texture];
            roughness *= evaluate_texture(roughness_texture, texturecoord)[1];
        }
        return roughness * roughness;
    }
}
vec3f evaluate_material_transmission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    auto transmission = material.transmission * make_shorter_vec(shape_color);
    if (material.transmission_texture >= 0) {
        auto& transmission_texture = scene.textures[material.transmission_texture];
        transmission *= make_shorter_vec(
            evaluate_texture(transmission_texture, texturecoord));
    }
    return transmission;
}
float evaluate_material_opacity(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    auto opacity = material.opacity * shape_color[3];
    if (material.opacity_texture >= 0) {
        auto& opacity_texture = scene.textures[material.opacity_texture];
        opacity *= evaluate_texture(opacity_texture, texturecoord)[3];
    }
    return opacity;
}

// Evaluates the microfacet_brdf at a location.
microfacet_brdf evaluate_material_brdf(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    auto brdf    = microfacet_brdf();
    brdf.diffuse = evaluate_material_diffuse(
        scene, material, texturecoord, shape_color);
    brdf.specular = evaluate_material_specular(
        scene, material, texturecoord, shape_color);
    brdf.transmission = evaluate_material_transmission(
        scene, material, texturecoord, shape_color);
    brdf.roughness = evaluate_material_roughness(
        scene, material, texturecoord, shape_color);
    brdf.refract = material.refract;
    if (brdf.diffuse != zero_vec3f) {
        brdf.roughness = clamp(brdf.roughness, 0.03f * 0.03f, 1.0f);
    } else if (brdf.roughness <= 0.03f * 0.03f)
        brdf.roughness = 0;
    return brdf;
}

bool is_brdf_delta(const microfacet_brdf& brdf) {
    return brdf.roughness == 0 && brdf.diffuse == zero_vec3f &&
           (brdf.specular != zero_vec3f || brdf.transmission != zero_vec3f);
}
bool is_brdf_zero(const microfacet_brdf& brdf) {
    return brdf.diffuse == zero_vec3f && brdf.specular == zero_vec3f &&
           brdf.transmission == zero_vec3f;
}

bool is_material_volume_homogeneus(const yocto_material& material) {
    return material.volume_density_texture < 0;
}
bool is_material_volume_colored(const yocto_material& material) {
    return !(material.volume_density[0] == material.volume_density[1] &&
             material.volume_density[1] == material.volume_density[2]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Merge scene into one another
void merge_scene(yocto_scene& merge_scene, const yocto_scene& merge_from) {
    log_error("this is  broken since we did not fix references");
    auto merge = [](auto& v1, auto& v2) {
        v1.insert(v1.end(), v2.begin(), v2.end());
    };
    auto merge_ = [](auto& v1, auto& v2) {
        v1.insert(v1.end(), v2.begin(), v2.end());
    };
    merge_(merge_scene.cameras, merge_from.cameras);
    merge(merge_scene.textures, merge_from.textures);
    merge(merge_scene.voltextures, merge_from.voltextures);
    merge(merge_scene.materials, merge_from.materials);
    merge(merge_scene.shapes, merge_from.shapes);
    merge(merge_scene.surfaces, merge_from.surfaces);
    merge(merge_scene.instances, merge_from.instances);
    merge(merge_scene.environments, merge_from.environments);
    merge_(merge_scene.nodes, merge_from.nodes);
    merge_(merge_scene.animations, merge_from.animations);
}

void print_stats(const yocto_scene& scene) {
    // using long long instead of uint64_t to avoid printf macros
    auto num_cameras      = (long long)0;
    auto num_shapes       = (long long)0;
    auto num_surfaces     = (long long)0;
    auto num_instances    = (long long)0;
    auto num_materials    = (long long)0;
    auto num_textures     = (long long)0;
    auto num_voltextures  = (long long)0;
    auto num_environments = (long long)0;
    auto num_nodes        = (long long)0;
    auto num_animations   = (long long)0;

    auto elem_points    = (long long)0;
    auto elem_lines     = (long long)0;
    auto elem_triangles = (long long)0;
    auto elem_quads     = (long long)0;
    auto vert_pos       = (long long)0;
    auto vert_norm      = (long long)0;
    auto vert_texcoord  = (long long)0;
    auto vert_color     = (long long)0;
    auto vert_radius    = (long long)0;
    auto vert_tangsp    = (long long)0;

    auto elem_quads_pos      = (long long)0;
    auto elem_quads_norm     = (long long)0;
    auto elem_quads_texcoord = (long long)0;
    auto vert_quads_pos      = (long long)0;
    auto vert_quads_norm     = (long long)0;
    auto vert_quads_texcoord = (long long)0;

    auto texel_hdr = (long long)0;
    auto texel_ldr = (long long)0;
    auto voxel_hdr = (long long)0;

    auto memory_imgs    = (long long)0;
    auto memory_vols    = (long long)0;
    auto memory_elems   = (long long)0;
    auto memory_verts   = (long long)0;
    auto memory_fvelems = (long long)0;
    auto memory_fvverts = (long long)0;

    auto bbox = compute_scene_bounds(scene);

    num_cameras      = scene.cameras.size();
    num_shapes       = scene.shapes.size();
    num_surfaces     = scene.surfaces.size();
    num_materials    = scene.materials.size();
    num_textures     = scene.textures.size();
    num_voltextures  = scene.voltextures.size();
    num_environments = scene.environments.size();
    num_instances    = scene.instances.size();
    num_nodes        = scene.nodes.size();
    num_animations   = scene.animations.size();

    for (auto& shape : scene.shapes) {
        elem_points += shape.points.size();
        elem_lines += shape.lines.size();
        elem_triangles += shape.triangles.size();
        elem_quads += shape.quads.size();
        vert_pos += shape.positions.size();
        vert_norm += shape.normals.size();
        vert_texcoord += shape.texturecoords.size();
        vert_color += shape.colors.size();
        vert_radius += shape.radius.size();
        vert_tangsp += shape.tangentspaces.size();
    }

    memory_elems = elem_points * sizeof(int) + elem_lines * sizeof(vec2i) +
                   elem_triangles * sizeof(vec3i) + elem_quads * sizeof(vec4i);
    memory_verts = vert_pos * sizeof(vec3f) + vert_norm * sizeof(vec3f) +
                   vert_texcoord * sizeof(vec2f) + vert_color * sizeof(vec4f) +
                   vert_tangsp * sizeof(vec4f) + vert_radius * sizeof(float);

    for (auto& surface : scene.surfaces) {
        elem_quads_pos += surface.quads_positions.size();
        elem_quads_norm += surface.quads_normals.size();
        elem_quads_texcoord += surface.quads_texturecoords.size();
        vert_quads_pos += surface.positions.size();
        vert_quads_norm += surface.normals.size();
        vert_quads_texcoord += surface.texturecoords.size();
    }

    memory_fvelems = elem_quads_pos * sizeof(vec4i) +
                     elem_quads_norm * sizeof(vec4i) +
                     elem_quads_texcoord * sizeof(vec4i);
    memory_fvverts = vert_quads_pos * sizeof(vec3f) +
                     vert_quads_norm * sizeof(vec3f) +
                     vert_quads_texcoord * sizeof(vec2f);

    for (auto& texture : scene.textures) {
        texel_hdr += texture.hdr_image.width() * texture.hdr_image.height();
        texel_ldr += texture.ldr_image.width() * texture.ldr_image.height();
    }
    memory_imgs = texel_hdr * sizeof(vec4f) + texel_ldr * sizeof(vec4b);

    for (auto& voltexture : scene.voltextures) {
        voxel_hdr += voltexture.volume_data.width() *
                     voltexture.volume_data.height() *
                     voltexture.volume_data.depth();
    }
    memory_vols = voxel_hdr * sizeof(float);

    printf("num_cameras: %lld\n", num_cameras);
    printf("num_shapes: %lld\n", num_shapes);
    printf("num_surface: %lld\n", num_surfaces);
    printf("num_instances: %lld\n", num_instances);
    printf("num_materials: %lld\n", num_materials);
    printf("num_textures: %lld\n", num_textures);
    printf("num_voltextures: %lld\n", num_voltextures);
    printf("num_environments: %lld\n", num_environments);
    printf("num_nodes: %lld\n", num_nodes);
    printf("num_animations: %lld\n", num_animations);

    printf("elem_points: %lld\n", elem_points);
    printf("elem_lines: %lld\n", elem_lines);
    printf("elem_triangles: %lld\n", elem_triangles);
    printf("elem_quads: %lld\n", elem_quads);
    printf("vert_pos: %lld\n", vert_pos);
    printf("vert_norm: %lld\n", vert_norm);
    printf("vert_texcoord: %lld\n", vert_texcoord);
    printf("vert_color: %lld\n", vert_color);
    printf("vert_radius: %lld\n", vert_radius);
    printf("vert_tangsp: %lld\n", vert_tangsp);

    printf("elem_points: %lld\n", elem_points);
    printf("elem_lines: %lld\n", elem_lines);
    printf("elem_triangles: %lld\n", elem_triangles);
    printf("elem_quads: %lld\n", elem_quads);
    printf("vert_pos: %lld\n", vert_pos);
    printf("vert_norm: %lld\n", vert_norm);
    printf("vert_texcoord: %lld\n", vert_texcoord);

    printf("texel_hdr: %lld\n", texel_hdr);
    printf("texel_ldr: %lld\n", texel_ldr);

    printf("memory_imgs: %lld\n", memory_imgs);
    printf("memory_vols: %lld\n", memory_vols);
    printf("memory_elems: %lld\n", memory_elems);
    printf("memory_verts: %lld\n", memory_verts);
    printf("memory_fvelems: %lld\n", memory_fvelems);
    printf("memory_fvverts: %lld\n", memory_fvverts);

    printf("bbox min: %g %g %g\n", bbox.min[0], bbox.min[1], bbox.min[2]);
    printf("bbox max: %g %g %g\n", bbox.max[0], bbox.max[1], bbox.max[2]);
}

}  // namespace yocto
