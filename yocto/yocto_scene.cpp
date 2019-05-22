//
// Implementation for Yocto/Scene.
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

#include "yocto_scene.h"
#include "yocto_random.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include <cassert>
#include <unordered_map>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {
using std::unordered_map;
}

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes a shape bounding box.
bbox3f compute_bounds(const yocto_shape& shape) {
    auto bbox = invalid_bbox3f;
    for (auto p : shape.positions) bbox += p;
    return bbox;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const yocto_scene& scene) {
    auto shape_bbox = vector<bbox3f>(scene.shapes.size());
    for (auto shape_id = 0; shape_id < scene.shapes.size(); shape_id++)
        shape_bbox[shape_id] = compute_bounds(scene.shapes[shape_id]);
    auto bbox = invalid_bbox3f;
    for (auto& instance : scene.instances) {
        bbox += transform_bbox(instance.frame, shape_bbox[instance.shape]);
    }
    return bbox;
}

// Compute vertex normals
void compute_normals(const yocto_shape& shape, vector<vec3f>& normals) {
    normals.assign(shape.positions.size(), {0, 0, 1});
    if (!shape.points.empty()) {
    } else if (!shape.lines.empty()) {
        compute_tangents(normals, shape.lines, shape.positions);
    } else if (!shape.triangles.empty()) {
        compute_normals(normals, shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
        compute_normals(normals, shape.quads, shape.positions);
    } else if (!shape.quads_positions.empty()) {
        compute_normals(normals, shape.quads_positions, shape.positions);
    } else {
        throw std::runtime_error("unknown element type");
    }
}

// Apply subdivision and displacement rules.
void subdivide_shape(yocto_shape& shape, int subdivision_level,
    bool catmull_clark, bool update_normals) {
    if (!subdivision_level) return;
    if (!shape.points.empty()) {
        throw runtime_error("point subdivision not supported");
    } else if (!shape.lines.empty()) {
        subdivide_lines(shape.lines, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius, subdivision_level);
    } else if (!shape.triangles.empty()) {
        subdivide_triangles(shape.triangles, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius, subdivision_level);
    } else if (!shape.quads.empty() && !catmull_clark) {
        subdivide_quads(shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius, subdivision_level);
    } else if (!shape.quads.empty() && catmull_clark) {
        subdivide_catmullclark(shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius, subdivision_level);
    } else if (!shape.quads_positions.empty() && !catmull_clark) {
        subdivide_quads(
            shape.quads_positions, shape.positions, subdivision_level);
        subdivide_quads(shape.quads_normals, shape.normals, subdivision_level);
        subdivide_quads(
            shape.quads_texcoords, shape.texcoords, subdivision_level);
    } else if (!shape.quads_positions.empty() && catmull_clark) {
        subdivide_catmullclark(
            shape.quads_positions, shape.positions, subdivision_level);
        subdivide_catmullclark(
            shape.quads_texcoords, shape.texcoords, subdivision_level, true);
    } else {
        throw runtime_error("empty shape");
    }

    if (update_normals) {
        if (!shape.quads_positions.empty()) {
            shape.quads_normals = shape.quads_positions;
        }
        compute_normals(shape, shape.normals);
    }
}
// Apply displacement to a shape
void displace_shape(yocto_shape& shape, const yocto_texture& displacement,
    float scale, bool update_normals) {
    if (shape.texcoords.empty()) {
        throw runtime_error("missing texture coordinates");
        return;
    }

    // simple case
    if (shape.quads_positions.empty()) {
        auto normals = shape.normals;
        if (shape.normals.empty()) compute_normals(shape, normals);
        for (auto vid = 0; vid < shape.positions.size(); vid++) {
            shape.positions[vid] +=
                normals[vid] * scale *
                mean(xyz(eval_texture(displacement, shape.texcoords[vid])));
        }
        if (update_normals || !shape.normals.empty()) {
            compute_normals(shape, shape.normals);
        }
    } else {
        // facevarying case
        auto offset = vector<float>(shape.positions.size(), 0);
        auto count  = vector<int>(shape.positions.size(), 0);
        for (auto fid = 0; fid < shape.quads_positions.size(); fid++) {
            auto qpos = shape.quads_positions[fid];
            auto qtxt = shape.quads_texcoords[fid];
            for (auto i = 0; i < 4; i++) {
                offset[qpos[i]] += scale * mean(xyz(eval_texture(displacement,
                                               shape.texcoords[qtxt[i]])));
                count[qpos[i]] += 1;
            }
        }
        auto normals = vector<vec3f>{shape.positions.size()};
        compute_normals(normals, shape.quads_positions, shape.positions);
        for (auto vid = 0; vid < shape.positions.size(); vid++) {
            shape.positions[vid] += normals[vid] * offset[vid] / count[vid];
        }
        if (update_normals || !shape.normals.empty()) {
            shape.quads_normals = shape.quads_positions;
            compute_normals(shape, shape.normals);
        }
    }
}

void tesselate_subdiv(yocto_scene& scene, yocto_subdiv& subdiv) {
    auto& shape           = scene.shapes[subdiv.tesselated_shape];
    shape.positions       = subdiv.positions;
    shape.normals         = subdiv.normals;
    shape.texcoords       = subdiv.texcoords;
    shape.colors          = subdiv.colors;
    shape.radius          = subdiv.radius;
    shape.points          = subdiv.points;
    shape.lines           = subdiv.lines;
    shape.triangles       = subdiv.triangles;
    shape.quads           = subdiv.quads;
    shape.quads_positions = subdiv.quads_positions;
    shape.quads_normals   = subdiv.quads_normals;
    shape.quads_texcoords = subdiv.quads_texcoords;
    shape.lines           = subdiv.lines;
    if (subdiv.subdivision_level) {
        subdivide_shape(shape, subdiv.subdivision_level, subdiv.catmull_clark,
            subdiv.compute_normals);
    }
    if (subdiv.displacement_texture >= 0) {
        displace_shape(shape, scene.textures[subdiv.displacement_texture],
            subdiv.displacement_scale, subdiv.compute_normals);
    }
}

// Updates tesselation.
void tesselate_subdivs(yocto_scene& scene) {
    for (auto& subdiv : scene.subdivs) tesselate_subdiv(scene, subdiv);
}

// Update animation transforms
void update_transforms(yocto_scene& scene, yocto_animation& animation,
    float time, const string& anim_group) {
    if (anim_group != "" && anim_group != animation.animation_group) return;

    if (!animation.translation_keyframes.empty()) {
        auto value = vec3f{0, 0, 0};
        switch (animation.interpolation_type) {
            case yocto_interpolation_type::step:
                value = keyframe_step(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = keyframe_linear(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = keyframe_bezier(animation.keyframes_times,
                    animation.translation_keyframes, time);
                break;
            default: throw runtime_error("should not have been here");
        }
        for (auto target : animation.node_targets)
            scene.nodes[target].translation = value;
    }
    if (!animation.rotation_keyframes.empty()) {
        auto value = vec4f{0, 0, 0, 1};
        switch (animation.interpolation_type) {
            case yocto_interpolation_type::step:
                value = keyframe_step(animation.keyframes_times,
                    animation.rotation_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = keyframe_linear(animation.keyframes_times,
                    animation.rotation_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = keyframe_bezier(animation.keyframes_times,
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
                value = keyframe_step(
                    animation.keyframes_times, animation.scale_keyframes, time);
                break;
            case yocto_interpolation_type::linear:
                value = keyframe_linear(
                    animation.keyframes_times, animation.scale_keyframes, time);
                break;
            case yocto_interpolation_type::bezier:
                value = keyframe_bezier(
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
    auto frame =
        parent * node.local * make_translation_frame(node.translation) *
        make_rotation_frame(node.rotation) * make_scaling_frame(node.scale);
    if (node.instance >= 0) scene.instances[node.instance].frame = frame;
    if (node.camera >= 0) scene.cameras[node.camera].frame = frame;
    if (node.environment >= 0)
        scene.environments[node.environment].frame = frame;
    for (auto child : node.children)
        update_transforms(scene, scene.nodes[child], frame);
}

// Update node transforms
void update_transforms(
    yocto_scene& scene, float time, const string& anim_group) {
    for (auto& agr : scene.animations)
        update_transforms(scene, agr, time, anim_group);
    for (auto& node : scene.nodes) node.children.clear();
    for (auto node_id = 0; node_id < scene.nodes.size(); node_id++) {
        auto& node = scene.nodes[node_id];
        if (node.parent >= 0)
            scene.nodes[node.parent].children.push_back(node_id);
    }
    for (auto& node : scene.nodes)
        if (node.parent < 0) update_transforms(scene, node);
}

// Compute animation range
vec2f compute_animation_range(
    const yocto_scene& scene, const string& anim_group) {
    if (scene.animations.empty()) return zero2f;
    auto range = vec2f{+float_max, -float_max};
    for (auto& animation : scene.animations) {
        if (anim_group != "" && animation.animation_group != anim_group)
            continue;
        range.x = min(range.x, animation.keyframes_times.front());
        range.y = max(range.y, animation.keyframes_times.back());
    }
    if (range.y < range.x) return zero2f;
    return range;
}

// Generate a distribution for sampling a shape uniformly based on area/length.
void sample_shape_cdf(const yocto_shape& shape, vector<float>& cdf) {
    cdf.clear();
    if (!shape.triangles.empty()) {
        sample_triangles_cdf(cdf, shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
        sample_quads_cdf(cdf, shape.quads, shape.positions);
    } else if (!shape.lines.empty()) {
        sample_lines_cdf(cdf, shape.lines, shape.positions);
    } else if (!shape.points.empty()) {
        sample_points_cdf(cdf, shape.points.size());
    } else if (!shape.quads_positions.empty()) {
        sample_quads_cdf(cdf, shape.quads_positions, shape.positions);
    } else {
        throw runtime_error("empty shape");
    }
}

// Sample a shape based on a distribution.
pair<int, vec2f> sample_shape(const yocto_shape& shape,
    const vector<float>& elements_cdf, float re, const vec2f& ruv) {
    // TODO: implement sampling without cdf
    if (elements_cdf.empty()) return {};
    if (!shape.triangles.empty()) {
        return sample_triangles(elements_cdf, re, ruv);
    } else if (!shape.quads.empty()) {
        return sample_quads(elements_cdf, re, ruv);
    } else if (!shape.lines.empty()) {
        return {sample_lines(elements_cdf, re, ruv.x).first, ruv};
    } else if (!shape.points.empty()) {
        return {sample_points(elements_cdf, re), ruv};
    } else if (!shape.quads_positions.empty()) {
        return sample_quads(elements_cdf, re, ruv);
    } else {
        return {0, zero2f};
    }
}

float sample_shape_pdf(const yocto_shape& shape,
    const vector<float>& elements_cdf, int element_id,
    const vec2f& element_uv) {
    // prob triangle * area triangle = area triangle mesh
    return 1 / elements_cdf.back();
}

// Update environment CDF for sampling.
void sample_environment_cdf(const yocto_scene& scene,
    const yocto_environment& environment, vector<float>& texels_cdf) {
    if (environment.emission_texture < 0) {
        texels_cdf.clear();
        return;
    }
    auto& texture = scene.textures[environment.emission_texture];
    auto  size    = texture_size(texture);
    texels_cdf.resize(size.x * size.y);
    if (size != zero2i) {
        for (auto i = 0; i < texels_cdf.size(); i++) {
            auto ij       = vec2i{i % size.x, i / size.x};
            auto th       = (ij.y + 0.5f) * pif / size.y;
            auto value    = lookup_texture(texture, ij.x, ij.y);
            texels_cdf[i] = max(xyz(value)) * sin(th);
            if (i) texels_cdf[i] += texels_cdf[i - 1];
        }
    } else {
        throw runtime_error("empty texture");
    }
}

// Sample an environment based on texels
vec3f sample_environment(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    float re, const vec2f& ruv) {
    if (!texels_cdf.empty() && environment.emission_texture >= 0) {
        auto& texture = scene.textures[environment.emission_texture];
        auto  idx     = sample_discrete(texels_cdf, re);
        auto  size    = texture_size(texture);
        auto  u       = (idx % size.x + 0.5f) / size.x;
        auto  v       = (idx / size.x + 0.5f) / size.y;
        return eval_direction(environment, {u, v});
    } else {
        return sample_sphere(ruv);
    }
}

// Sample an environment based on texels
float sample_environment_pdf(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    const vec3f& direction) {
    if (!texels_cdf.empty() && environment.emission_texture >= 0) {
        auto& texture  = scene.textures[environment.emission_texture];
        auto  size     = texture_size(texture);
        auto  texcoord = eval_texcoord(environment, direction);
        auto  i        = (int)(texcoord.x * size.x);
        auto  j        = (int)(texcoord.y * size.y);
        auto  idx      = j * size.x + i;
        auto  prob  = sample_discrete_pdf(texels_cdf, idx) / texels_cdf.back();
        auto  angle = (2 * pif / size.x) * (pif / size.y) *
                     sin(pif * (j + 0.5f) / size.y);
        return prob / angle;
    } else {
        return sample_sphere_pdf(direction);
    }
}

// Add missing names and resolve duplicated names.
void normalize_uris(yocto_scene& scene) {
    auto normalize = [](string& name, const string& base, const string& ext,
                         int num) {
        for (auto& c : name) {
            if (c == ':' || c == ' ') c = '_';
        }
        if (name.empty()) name = base + "_" + to_string(num);
        if (get_dirname(name).empty()) name = base + "s/" + name;
        if (get_extension(name).empty()) name = name + "." + ext;
    };
    for (auto id = 0; id < scene.cameras.size(); id++)
        normalize(scene.cameras[id].uri, "camera", "yaml", id);
    for (auto id = 0; id < scene.textures.size(); id++)
        normalize(scene.textures[id].uri, "texture", "png", id);
    for (auto id = 0; id < scene.voltextures.size(); id++)
        normalize(scene.voltextures[id].uri, "volume", "yvol", id);
    for (auto id = 0; id < scene.materials.size(); id++)
        normalize(scene.materials[id].uri, "material", "yaml", id);
    for (auto id = 0; id < scene.shapes.size(); id++)
        normalize(scene.shapes[id].uri, "shape", "ply", id);
    for (auto id = 0; id < scene.instances.size(); id++)
        normalize(scene.instances[id].uri, "instance", "yaml", id);
    for (auto id = 0; id < scene.animations.size(); id++)
        normalize(scene.animations[id].uri, "animation", "yaml", id);
    for (auto id = 0; id < scene.nodes.size(); id++)
        normalize(scene.nodes[id].uri, "node", "yaml", id);
}
void rename_instances(yocto_scene& scene) {
    auto shape_names = vector<string>(scene.shapes.size());
    for (auto sid = 0; sid < scene.shapes.size(); sid++) {
        shape_names[sid] = get_basename(scene.shapes[sid].uri);
    }
    auto shape_count = vector<vec2i>(scene.shapes.size(), vec2i{0, 0});
    for (auto& instance : scene.instances) shape_count[instance.shape].y += 1;
    for (auto& instance : scene.instances) {
        if (shape_count[instance.shape].y == 1) {
            instance.uri = format(
                "instances/{}.yaml", shape_names[instance.shape]);
        } else {
            instance.uri = format("instances/{}-{:{}}.yaml",
                shape_names[instance.shape], shape_count[instance.shape].x++,
                (int)ceil(log10(shape_count[instance.shape].y)));
        }
    }
}

// Normalized a scaled color in a material
void normalize_scaled_color(float& scale, vec3f& color) {
    auto scaled = scale * color;
    if (max(scaled) == 0) {
        scale = 0;
        color = {1, 1, 1};
    } else {
        scale = max(scaled);
        color = scaled / max(scaled);
    }
}

// Add missing tangent space if needed.
void add_tangent_spaces(yocto_scene& scene) {
    for (auto& instance : scene.instances) {
        auto& material = scene.materials[instance.material];
        if (material.normal_texture < 0) continue;
        auto& shape = scene.shapes[instance.shape];
        if (!shape.tangents.empty() || shape.texcoords.empty()) continue;
        if (!shape.triangles.empty()) {
            if (shape.normals.empty()) {
                shape.normals.resize(shape.positions.size());
                compute_normals(
                    shape.normals, shape.triangles, shape.positions);
            }
            shape.tangents.resize(shape.positions.size());
            compute_tangent_spaces(shape.tangents, shape.triangles,
                shape.positions, shape.normals, shape.texcoords);
        } else {
            throw runtime_error("type not supported");
        }
    }
}

// Add missing materials.
void add_materials(yocto_scene& scene) {
    auto material_id = -1;
    for (auto& instance : scene.instances) {
        if (instance.material >= 0) continue;
        if (material_id < 0) {
            auto material    = yocto_material{};
            material.uri     = "materails/default.yaml";
            material.diffuse = {0.2f, 0.2f, 0.2f};
            scene.materials.push_back(material);
            material_id = (int)scene.materials.size() - 1;
        }
        instance.material = material_id;
    }
}

// Add missing radius.
void add_radius(yocto_scene& scene, float radius) {
    for (auto& shape : scene.shapes) {
        if (shape.points.empty() && shape.lines.empty()) continue;
        if (!shape.radius.empty()) continue;
        shape.radius.assign(shape.positions.size(), radius);
    }
}

// Add missing cameras.
void add_cameras(yocto_scene& scene) {
    if (scene.cameras.empty()) {
        auto camera = yocto_camera{};
        camera.uri  = "cameras/default.yaml";
        set_view(camera, compute_bounds(scene), {0, 0, 1});
        scene.cameras.push_back(camera);
    }
}

// Add a sky environment
void add_sky(yocto_scene& scene, float sun_angle) {
    auto texture = yocto_texture{};
    texture.uri  = "textures/sky.hdr";
    make_sunsky(texture.hdr_image, {1024, 512}, sun_angle);
    scene.textures.push_back(texture);
    auto environment             = yocto_environment{};
    environment.uri              = "environments/default.yaml";
    environment.emission         = {1, 1, 1};
    environment.emission_texture = (int)scene.textures.size() - 1;
    scene.environments.push_back(environment);
}

// Reduce memory usage
void trim_memory(yocto_scene& scene) {
    for (auto& shape : scene.shapes) {
        shape.points.shrink_to_fit();
        shape.lines.shrink_to_fit();
        shape.triangles.shrink_to_fit();
        shape.quads.shrink_to_fit();
        shape.quads_positions.shrink_to_fit();
        shape.quads_normals.shrink_to_fit();
        shape.quads_texcoords.shrink_to_fit();
        shape.positions.shrink_to_fit();
        shape.normals.shrink_to_fit();
        shape.texcoords.shrink_to_fit();
        shape.colors.shrink_to_fit();
        shape.radius.shrink_to_fit();
        shape.tangents.shrink_to_fit();
    }
    for (auto& texture : scene.textures) {
        texture.ldr_image.shrink_to_fit();
        texture.hdr_image.shrink_to_fit();
    }
    scene.cameras.shrink_to_fit();
    scene.shapes.shrink_to_fit();
    scene.instances.shrink_to_fit();
    scene.materials.shrink_to_fit();
    scene.textures.shrink_to_fit();
    scene.environments.shrink_to_fit();
    scene.voltextures.shrink_to_fit();
    scene.nodes.shrink_to_fit();
    scene.animations.shrink_to_fit();
}

// Checks for validity of the scene.
vector<string> validate_scene(const yocto_scene& scene, bool skip_textures) {
    auto errs        = vector<string>();
    auto check_names = [&errs](const auto& vals, const string& base) {
        auto used = unordered_map<string, int>();
        used.reserve(vals.size());
        for (auto& value : vals) used[value.uri] += 1;
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
                errs.push_back("empty texture " + value.uri);
            }
        }
    };

    check_names(scene.cameras, "camera");
    check_names(scene.shapes, "shape");
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
void print_validation(const yocto_scene& scene, bool skip_textures) {
    for (auto err : validate_scene(scene, skip_textures))
        printf("%s [validation]\n", err.c_str());
}

void build_bvh(
    bvh_scene& bvh, const yocto_scene& scene, const bvh_params& params) {
    bvh.shapes.resize(scene.shapes.size());
    for (auto idx = 0; idx < scene.shapes.size(); idx++) {
        auto& shape = scene.shapes[idx];
        auto& sbvh  = bvh.shapes[idx];
#if YOCTO_EMBREE
        // call Embree if needed
        if (params.use_embree) {
            if (params.embree_compact &&
                shape.positions.size() == shape.positions.capacity()) {
                ((yocto_shape&)shape)
                    .positions.reserve(shape.positions.size() + 1);
            }
        }
#endif
        sbvh.points          = shape.points;
        sbvh.lines           = shape.lines;
        sbvh.triangles       = shape.triangles;
        sbvh.quads           = shape.quads;
        sbvh.quads_positions = shape.quads_positions;
        sbvh.positions       = shape.positions;
        sbvh.radius          = shape.radius;
    }
    if (!scene.instances.empty()) {
        bvh.instances = {&scene.instances[0].frame, (int)scene.instances.size(),
            sizeof(scene.instances[0])};
    } else {
        bvh.instances = {};
    }

    build_bvh(bvh, params);
}

void refit_bvh(bvh_scene& bvh, const yocto_scene& scene,
    const vector<int>& updated_shapes, const bvh_params& params) {
    for (auto idx : updated_shapes) {
        auto& shape = scene.shapes[idx];
        auto& sbvh  = bvh.shapes[idx];
#if YOCTO_EMBREE
        // call Embree if needed
        if (params.use_embree) {
            if (params.embree_compact &&
                shape.positions.size() == shape.positions.capacity()) {
                ((yocto_shape&)shape)
                    .positions.reserve(shape.positions.size() + 1);
            }
        }
#endif
        sbvh.points          = shape.points;
        sbvh.lines           = shape.lines;
        sbvh.triangles       = shape.triangles;
        sbvh.quads           = shape.quads;
        sbvh.quads_positions = shape.quads_positions;
        sbvh.positions       = shape.positions;
        sbvh.radius          = shape.radius;
    }
    if (!scene.instances.empty()) {
        bvh.instances = {&scene.instances[0].frame, (int)scene.instances.size(),
            sizeof(scene.instances[0])};
    } else {
        bvh.instances = {};
    }

    refit_bvh(bvh, updated_shapes, params);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Shape element normal.
vec3f eval_element_normal(const yocto_shape& shape, int element_id) {
    auto norm = zero3f;
    if (!shape.triangles.empty()) {
        auto t = shape.triangles[element_id];
        norm   = triangle_normal(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
    } else if (!shape.quads.empty()) {
        auto q = shape.quads[element_id];
        norm   = quad_normal(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
    } else if (!shape.lines.empty()) {
        auto l = shape.lines[element_id];
        norm   = line_tangent(shape.positions[l.x], shape.positions[l.y]);
    } else if (!shape.quads_positions.empty()) {
        auto q = shape.quads_positions[element_id];
        norm   = quad_normal(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
    } else {
        throw runtime_error("empty shape");
        norm = {0, 0, 1};
    }
    return norm;
}

// Shape element normal.
pair<vec3f, bool> eval_element_tangents(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (!shape.triangles.empty()) {
        auto t    = shape.triangles[element_id];
        auto norm = triangle_normal(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
        auto txty = pair<vec3f, vec3f>();
        if (shape.texcoords.empty()) {
            txty = triangle_tangents_fromuv(shape.positions[t.x],
                shape.positions[t.y], shape.positions[t.z], {0, 0}, {1, 0},
                {0, 1});
        } else {
            txty = triangle_tangents_fromuv(shape.positions[t.x],
                shape.positions[t.y], shape.positions[t.z],
                shape.texcoords[t.x], shape.texcoords[t.y],
                shape.texcoords[t.z]);
        }
        auto tx = txty.first, ty = txty.second;
        tx     = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        return {tx, s};
    } else if (!shape.quads.empty()) {
        auto q    = shape.quads[element_id];
        auto norm = quad_normal(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
        auto txty = pair<vec3f, vec3f>();
        if (shape.texcoords.empty()) {
            txty = quad_tangents_fromuv(shape.positions[q.x],
                shape.positions[q.y], shape.positions[q.z],
                shape.positions[q.w], {0, 0}, {1, 0}, {0, 1}, {1, 1},
                element_uv);
        } else {
            txty = quad_tangents_fromuv(shape.positions[q.x],
                shape.positions[q.y], shape.positions[q.z],
                shape.positions[q.w], shape.texcoords[q.x],
                shape.texcoords[q.y], shape.texcoords[q.z],
                shape.texcoords[q.w], element_uv);
        }
        auto tx = txty.first, ty = txty.second;
        tx     = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        return {tx, s};
    } else if (!shape.quads_positions.empty()) {
        auto q    = shape.quads_positions[element_id];
        auto norm = quad_normal(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
        auto txty = pair<vec3f, vec3f>();
        if (shape.texcoords.empty()) {
            txty = quad_tangents_fromuv(shape.positions[q.x],
                shape.positions[q.y], shape.positions[q.z],
                shape.positions[q.w], {0, 0}, {1, 0}, {0, 1}, {1, 1},
                element_uv);
        } else {
            auto qt = shape.quads_texcoords[element_id];
            txty    = quad_tangents_fromuv(shape.positions[q.x],
                shape.positions[q.y], shape.positions[q.z],
                shape.positions[q.w], shape.texcoords[qt.x],
                shape.texcoords[qt.y], shape.texcoords[qt.z],
                shape.texcoords[qt.w], element_uv);
        }
        auto tx = txty.first, ty = txty.second;
        tx     = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        return {tx, s};
    } else {
        return {zero3f, false};
    }
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T evaluate_shape_elem(const yocto_shape& shape,
    const vector<vec4i>& facevarying_quads, const vector<T>& vals,
    int element_id, const vec2f& element_uv) {
    if (vals.empty()) return {};
    if (!shape.triangles.empty()) {
        auto t = shape.triangles[element_id];
        return interpolate_triangle(
            vals[t.x], vals[t.y], vals[t.z], element_uv);
    } else if (!shape.quads.empty()) {
        auto q = shape.quads[element_id];
        if (q.w == q.z)
            return interpolate_triangle(
                vals[q.x], vals[q.y], vals[q.z], element_uv);
        return interpolate_quad(
            vals[q.x], vals[q.y], vals[q.z], vals[q.w], element_uv);
    } else if (!shape.lines.empty()) {
        auto l = shape.lines[element_id];
        return interpolate_line(vals[l.x], vals[l.y], element_uv.x);
    } else if (!shape.points.empty()) {
        return vals[shape.points[element_id]];
    } else if (!shape.quads_positions.empty()) {
        auto q = facevarying_quads[element_id];
        if (q.w == q.z)
            return interpolate_triangle(
                vals[q.x], vals[q.y], vals[q.z], element_uv);
        return interpolate_quad(
            vals[q.x], vals[q.y], vals[q.z], vals[q.w], element_uv);
    } else {
        return {};
    }
}

// Shape values interpolated using barycentric coordinates
vec3f eval_position(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    return evaluate_shape_elem(
        shape, shape.quads_positions, shape.positions, element_id, element_uv);
}
vec3f eval_normal(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.normals.empty()) return eval_element_normal(shape, element_id);
    return normalize(evaluate_shape_elem(
        shape, shape.quads_normals, shape.normals, element_id, element_uv));
}
vec2f eval_texcoord(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.texcoords.empty()) return element_uv;
    return evaluate_shape_elem(
        shape, shape.quads_texcoords, shape.texcoords, element_id, element_uv);
}
vec4f eval_color(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.colors.empty()) return {1, 1, 1, 1};
    return evaluate_shape_elem(shape, {}, shape.colors, element_id, element_uv);
}
float eval_radius(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.radius.empty()) return 0.001f;
    return evaluate_shape_elem(shape, {}, shape.radius, element_id, element_uv);
}
pair<vec3f, bool> eval_tangsp(
    const yocto_shape& shape, int element_id, const vec2f& element_uv) {
    if (shape.tangents.empty())
        return eval_element_tangents(shape, element_id, element_uv);
    auto tangsp = evaluate_shape_elem(
        shape, {}, shape.tangents, element_id, element_uv);
    return {xyz(tangsp), tangsp.w < 0};
}
// Shading normals including material perturbations.
vec3f eval_perturbed_normal(const yocto_scene& scene, const yocto_shape& shape,
    int element_id, const vec2f& element_uv, const vec3f& normalmap) {
    auto normal = eval_normal(shape, element_id, element_uv);
    if (shape.triangles.empty() && shape.quads.empty()) return normal;
    auto [tu, left_handed] = eval_tangsp(shape, element_id, element_uv);
    tu                     = orthonormalize(tu, normal);
    auto tv = normalize(cross(normal, tu) * (left_handed ? -1.0f : 1.0f));
    normal  = normalize(
        normalmap.x * tu + normalmap.y * tv + normalmap.z * normal);
    return normal;
}

// Instance values interpolated using barycentric coordinates.
vec3f eval_position(const yocto_scene& scene, const yocto_instance& instance,
    int element_id, const vec2f& element_uv) {
    return transform_point(instance.frame,
        eval_position(scene.shapes[instance.shape], element_id, element_uv));
}
vec3f eval_normal(const yocto_scene& scene, const yocto_instance& instance,
    int element_id, const vec2f& element_uv, bool non_rigid_frame) {
    auto normal = eval_normal(
        scene.shapes[instance.shape], element_id, element_uv);
    return non_rigid_frame
               ? transform_normal((const affine3f&)instance.frame, normal)
               : transform_normal(instance.frame, normal);
}
vec3f eval_shading_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    const vec3f& direction, bool non_rigid_frame) {
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    if (!shape.points.empty()) {
        return -direction;
    } else if (!shape.lines.empty()) {
        auto normal = eval_normal(
            scene, instance, element_id, element_uv, non_rigid_frame);
        return orthonormalize(-direction, normal);
    } else if (material.normal_texture < 0) {
        return eval_normal(
            scene, instance, element_id, element_uv, non_rigid_frame);
    } else {
        auto& normal_texture = scene.textures[material.normal_texture];
        auto  normalmap =
            xyz(eval_texture(normal_texture,
                eval_texcoord(shape, element_id, element_uv), true)) *
                2 -
            1;
        normalmap.y = -normalmap.y;  // flip vertical axis
        auto normal = eval_perturbed_normal(scene, scene.shapes[instance.shape],
            element_id, element_uv, normalmap);
        return non_rigid_frame
                   ? transform_normal((const affine3f&)instance.frame, normal)
                   : transform_normal(instance.frame, normal);
    }
}
// Instance element values.
vec3f eval_element_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, bool non_rigid_frame) {
    auto normal = eval_element_normal(scene.shapes[instance.shape], element_id);
    return non_rigid_frame
               ? transform_normal((const affine3f&)instance.frame, normal)
               : transform_normal(instance.frame, normal);
}
// Instance material
material_point eval_material(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv) {
    auto& shape     = scene.shapes[instance.shape];
    auto& material  = scene.materials[instance.material];
    auto  texcoords = eval_texcoord(shape, element_id, element_uv);
    auto  color     = eval_color(shape, element_id, element_uv);
    return eval_material(scene, material, texcoords, color);
}

// Environment texture coordinates from the direction.
vec2f eval_texcoord(
    const yocto_environment& environment, const vec3f& direction) {
    auto wl = transform_direction_inverse(environment.frame, direction);
    auto environment_uv = vec2f{
        atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
    if (environment_uv.x < 0) environment_uv.x += 1;
    return environment_uv;
}
// Evaluate the environment direction.
vec3f eval_direction(
    const yocto_environment& environment, const vec2f& environment_uv) {
    return transform_direction(environment.frame,
        {cos(environment_uv.x * 2 * pif) * sin(environment_uv.y * pif),
            cos(environment_uv.y * pif),
            sin(environment_uv.x * 2 * pif) * sin(environment_uv.y * pif)});
}
// Evaluate the environment color.
vec3f eval_environment(const yocto_scene& scene,
    const yocto_environment& environment, const vec3f& direction) {
    auto emission = environment.emission;
    if (environment.emission_texture >= 0) {
        auto& emission_texture = scene.textures[environment.emission_texture];
        emission *= xyz(eval_texture(
            emission_texture, eval_texcoord(environment, direction)));
    }
    return emission;
}
// Evaluate all environment color.
vec3f eval_environment(const yocto_scene& scene, const vec3f& direction) {
    auto emission = zero3f;
    for (auto& environment : scene.environments)
        emission += eval_environment(scene, environment, direction);
    return emission;
}

// Check texture size
vec2i texture_size(const yocto_texture& texture) {
    if (!texture.hdr_image.empty()) {
        return texture.hdr_image.size();
    } else if (!texture.ldr_image.empty()) {
        return texture.ldr_image.size();
    } else {
        return zero2i;
    }
}

// Lookup a texture value
vec4f lookup_texture(
    const yocto_texture& texture, int i, int j, bool ldr_as_linear) {
    if (!texture.hdr_image.empty()) {
        return texture.hdr_image[{i, j}];
    } else if (!texture.ldr_image.empty() && !ldr_as_linear) {
        return srgb_to_rgb(byte_to_float(texture.ldr_image[{i, j}]));
    } else if (!texture.ldr_image.empty() && ldr_as_linear) {
        return byte_to_float(texture.ldr_image[{i, j}]);
    } else {
        return zero4f;
    }
}

// Evaluate a texture
vec4f eval_texture(const yocto_texture& texture, const vec2f& texcoord,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
    if (texture.hdr_image.empty() && texture.ldr_image.empty())
        return {1, 1, 1, 1};

    // get image width/height
    auto size  = texture_size(texture);
    auto width = size.x, height = size.y;

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (clamp_to_edge) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * width;
        t = clamp(texcoord.y, 0.0f, 1.0f) * height;
    } else {
        s = fmod(texcoord.x, 1.0f) * width;
        if (s < 0) s += width;
        t = fmod(texcoord.y, 1.0f) * height;
        if (t < 0) t += height;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, width - 1), j = clamp((int)t, 0, height - 1);
    auto ii = (i + 1) % width, jj = (j + 1) % height;
    auto u = s - i, v = t - j;

    if (no_interpolation) return lookup_texture(texture, i, j, ldr_as_linear);

    // handle interpolation
    return lookup_texture(texture, i, j, ldr_as_linear) * (1 - u) * (1 - v) +
           lookup_texture(texture, i, jj, ldr_as_linear) * (1 - u) * v +
           lookup_texture(texture, ii, j, ldr_as_linear) * u * (1 - v) +
           lookup_texture(texture, ii, jj, ldr_as_linear) * u * v;
}

// Lookup a texture value
float lookup_voltexture(
    const yocto_voltexture& texture, int i, int j, int k, bool ldr_as_linear) {
    if (!texture.volume_data.empty()) {
        return texture.volume_data[{i, j, k}];
    } else {
        return 0;
    }
}

// Evaluate a volume texture
float eval_voltexture(const yocto_voltexture& texture, const vec3f& texcoord,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
    if (texture.volume_data.empty()) return 1;

    // get image width/height
    auto width  = texture.volume_data.size().x;
    auto height = texture.volume_data.size().y;
    auto depth  = texture.volume_data.size().z;

    // get coordinates normalized for tiling
    auto s = clamp((texcoord.x + 1.0f) * 0.5f, 0.0f, 1.0f) * width;
    auto t = clamp((texcoord.y + 1.0f) * 0.5f, 0.0f, 1.0f) * height;
    auto r = clamp((texcoord.z + 1.0f) * 0.5f, 0.0f, 1.0f) * depth;

    // get image coordinates and residuals
    auto i  = clamp((int)s, 0, width - 1);
    auto j  = clamp((int)t, 0, height - 1);
    auto k  = clamp((int)r, 0, depth - 1);
    auto ii = (i + 1) % width, jj = (j + 1) % height, kk = (k + 1) % depth;
    auto u = s - i, v = t - j, w = r - k;

    // nearest-neighbor interpolation
    if (no_interpolation) {
        i = u < 0.5 ? i : min(i + 1, width - 1);
        j = v < 0.5 ? j : min(j + 1, height - 1);
        k = w < 0.5 ? k : min(k + 1, depth - 1);
        return lookup_voltexture(texture, i, j, k, ldr_as_linear);
    }

    // trilinear interpolation
    return lookup_voltexture(texture, i, j, k, ldr_as_linear) * (1 - u) *
               (1 - v) * (1 - w) +
           lookup_voltexture(texture, ii, j, k, ldr_as_linear) * u * (1 - v) *
               (1 - w) +
           lookup_voltexture(texture, i, jj, k, ldr_as_linear) * (1 - u) * v *
               (1 - w) +
           lookup_voltexture(texture, i, j, kk, ldr_as_linear) * (1 - u) *
               (1 - v) * w +
           lookup_voltexture(texture, i, jj, kk, ldr_as_linear) * (1 - u) * v *
               w +
           lookup_voltexture(texture, ii, j, kk, ldr_as_linear) * u * (1 - v) *
               w +
           lookup_voltexture(texture, ii, jj, k, ldr_as_linear) * u * v *
               (1 - w) +
           lookup_voltexture(texture, ii, jj, kk, ldr_as_linear) * u * v * w;
}

// Set and evaluate camera parameters. Setters take zeros as default values.
float camera_fovx(const yocto_camera& camera) {
    assert(!camera.orthographic);
    return 2 * atan(camera.film_width / (2 * camera.focal_length));
}
float camera_fovy(const yocto_camera& camera) {
    assert(!camera.orthographic);
    return 2 * atan(camera.film_height / (2 * camera.focal_length));
}
float camera_aspect(const yocto_camera& camera) {
    return camera.film_width / camera.film_height;
}
vec2i camera_image_size(const yocto_camera& camera, const vec2i& size_) {
    auto size = size_;
    if (size == zero2i) size = {1280, 720};
    if (size.x != 0 && size.y != 0) {
        if (size.x * camera.film_height / camera.film_width > size.y) {
            size.x = 0;
        } else {
            size.y = 0;
        }
    }
    if (size.x == 0) {
        size.x = (int)round(size.y * camera.film_width / camera.film_height);
    }
    if (size.y == 0) {
        size.y = (int)round(size.x * camera.film_height / camera.film_width);
    }
    return size;
}
void set_perspectivey(
    yocto_camera& camera, float fovy, float aspect, float focus, float height) {
    camera.orthographic   = false;
    camera.film_width     = height * aspect;
    camera.film_height    = height;
    camera.focus_distance = focus;
    auto distance         = camera.film_height / (2 * tan(fovy / 2));
    if (focus < float_max) {
        camera.focal_length = camera.focus_distance * distance /
                              (camera.focus_distance + distance);
    } else {
        camera.focal_length = distance;
    }
}
void set_perspectivex(
    yocto_camera& camera, float fovx, float aspect, float focus, float width) {
    camera.orthographic   = false;
    camera.film_width     = width;
    camera.film_height    = width / aspect;
    camera.focus_distance = focus;
    auto distance         = camera.film_width / (2 * tan(fovx / 2));
    if (focus < float_max) {
        camera.focal_length = camera.focus_distance * distance /
                              (camera.focus_distance + distance);
    } else {
        camera.focal_length = distance;
    }
}

// add missing camera
void set_view(yocto_camera& camera, const bbox3f& bbox,
    const vec3f& view_direction, float width, float height, float focal) {
    camera.orthographic = false;
    if (width != 0) camera.film_width = width;
    if (height != 0) camera.film_height = height;
    if (focal != 0) camera.focal_length = focal;
    auto bbox_center = (bbox.max + bbox.min) / 2.0f;
    auto bbox_radius = length(bbox.max - bbox.min) / 2;
    auto camera_dir  = (view_direction == zero3f) ? camera.frame.o - bbox_center
                                                 : view_direction;
    if (camera_dir == zero3f) camera_dir = {0, 0, 1};
    auto camera_fov = min(camera_fovx(camera), camera_fovy(camera));
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
ray3f eval_perspective_camera(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
    auto distance = camera.focal_length;
    if (camera.focus_distance < float_max) {
        distance = camera.focal_length * camera.focus_distance /
                   (camera.focus_distance - camera.focal_length);
    }
    if (camera.lens_aperture) {
        auto e = vec3f{(lens_uv.x - 0.5f) * camera.lens_aperture,
            (lens_uv.y - 0.5f) * camera.lens_aperture, 0};
        auto q = vec3f{camera.film_width * (0.5f - image_uv.x),
            camera.film_height * (image_uv.y - 0.5f), distance};
        // distance of the image of the point
        auto distance1 = camera.focal_length * distance /
                         (distance - camera.focal_length);
        auto q1 = -q * distance1 / distance;
        auto d  = normalize(q1 - e);
        // auto q1 = - normalize(q) * camera.focus_distance / normalize(q).z;
        auto ray = make_ray(transform_point(camera.frame, e),
            transform_direction(camera.frame, d));
        return ray;
    } else {
        auto e   = zero3f;
        auto q   = vec3f{camera.film_width * (0.5f - image_uv.x),
            camera.film_height * (image_uv.y - 0.5f), distance};
        auto q1  = -q;
        auto d   = normalize(q1 - e);
        auto ray = make_ray(transform_point(camera.frame, e),
            transform_direction(camera.frame, d));
        return ray;
    }
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_orthographic_camera(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
    if (camera.lens_aperture) {
        auto scale = 1 / camera.focal_length;
        auto q     = vec3f{camera.film_width * (0.5f - image_uv.x) * scale,
            camera.film_height * (image_uv.y - 0.5f) * scale, scale};
        auto q1    = vec3f{-q.x, -q.y, -camera.focus_distance};
        auto e     = vec3f{-q.x, -q.y, 0} +
                 vec3f{(lens_uv.x - 0.5f) * camera.lens_aperture,
                     (lens_uv.y - 0.5f) * camera.lens_aperture, 0};
        auto d   = normalize(q1 - e);
        auto ray = make_ray(transform_point(camera.frame, e),
            transform_direction(camera.frame, d));
        return ray;
    } else {
        auto scale = 1 / camera.focal_length;
        auto q     = vec3f{camera.film_width * (0.5f - image_uv.x) * scale,
            camera.film_height * (image_uv.y - 0.5f) * scale, scale};
        auto q1    = -q;
        auto e     = vec3f{-q.x, -q.y, 0};
        auto d     = normalize(q1 - e);
        auto ray   = make_ray(transform_point(camera.frame, e),
            transform_direction(camera.frame, d));
        return ray;
    }
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
    if (camera.orthographic)
        return eval_orthographic_camera(camera, image_uv, lens_uv);
    else
        return eval_perspective_camera(camera, image_uv, lens_uv);
}

// Generates a ray from a camera.
ray3f eval_camera(const yocto_camera& camera, const vec2i& image_ij,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv) {
    auto image_uv = vec2f{(image_ij.x + pixel_uv.x) / image_size.x,
        (image_ij.y + pixel_uv.y) / image_size.y};
    return eval_camera(camera, image_uv, lens_uv);
}

// Generates a ray from a camera.
ray3f eval_camera(const yocto_camera& camera, int idx, const vec2i& image_size,
    const vec2f& pixel_uv, const vec2f& lens_uv) {
    auto image_ij = vec2i{idx % image_size.x, idx / image_size.x};
    auto image_uv = vec2f{(image_ij.x + pixel_uv.x) / image_size.x,
        (image_ij.y + pixel_uv.y) / image_size.y};
    return eval_camera(camera, image_uv, lens_uv);
}

// Evaluates the microfacet_brdf at a location.
material_point eval_material(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color) {
    auto point = material_point{};
    // factors
    point.emission       = material.emission * xyz(shape_color);
    point.diffuse        = material.diffuse * xyz(shape_color);
    point.specular       = material.specular;
    auto metallic        = material.metallic;
    point.roughness      = material.roughness;
    point.coat           = material.coat;
    point.transmission   = material.transmission;
    auto voltransmission = material.voltransmission;
    point.volemission    = material.volemission;
    point.volscatter     = material.volscatter;
    point.volanisotropy  = material.volanisotropy;
    auto volscale        = material.volscale;
    point.opacity        = material.opacity * shape_color.w;
    point.thin           = material.thin;

    // textures
    if (material.emission_texture >= 0) {
        auto& emission_texture = scene.textures[material.emission_texture];
        point.emission *= xyz(eval_texture(emission_texture, texturecoord));
    }
    if (material.diffuse_texture >= 0) {
        auto& diffuse_texture = scene.textures[material.diffuse_texture];
        auto  base_txt        = eval_texture(diffuse_texture, texturecoord);
        point.diffuse *= xyz(base_txt);
        point.opacity *= base_txt.w;
    }
    if (material.metallic_texture >= 0) {
        auto& metallic_texture = scene.textures[material.metallic_texture];
        auto  metallic_txt     = eval_texture(metallic_texture, texturecoord);
        metallic *= metallic_txt.z;
        if (material.gltf_textures) {
            point.roughness *= metallic_txt.x;
        }
    }
    if (material.specular_texture >= 0) {
        auto& specular_texture = scene.textures[material.specular_texture];
        auto  specular_txt     = eval_texture(specular_texture, texturecoord);
        point.specular *= xyz(specular_txt);
        if (material.gltf_textures) {
            auto glossiness = 1 - point.roughness;
            glossiness *= specular_txt.w;
            point.roughness = 1 - glossiness;
        }
    }
    if (material.roughness_texture >= 0) {
        auto& roughness_texture = scene.textures[material.roughness_texture];
        point.roughness *= eval_texture(roughness_texture, texturecoord).x;
    }
    if (material.transmission_texture >= 0) {
        auto& transmission_texture =
            scene.textures[material.transmission_texture];
        point.transmission *= xyz(
            eval_texture(transmission_texture, texturecoord));
    }
    if (material.subsurface_texture >= 0) {
        auto& subsurface_texture = scene.textures[material.subsurface_texture];
        point.volscatter *= xyz(eval_texture(subsurface_texture, texturecoord));
    }
    if (material.opacity_texture >= 0) {
        auto& opacity_texture = scene.textures[material.opacity_texture];
        point.opacity *= mean(xyz(eval_texture(opacity_texture, texturecoord)));
    }
    if (material.coat_texture >= 0) {
        auto& coat_texture = scene.textures[material.coat_texture];
        point.coat *= xyz(eval_texture(coat_texture, texturecoord));
    }
    if (metallic) {
        point.specular = point.specular * (1 - metallic) +
                         metallic * point.diffuse;
        point.diffuse = metallic * point.diffuse * (1 - metallic);
    }
    if (point.diffuse != zero3f || point.roughness) {
        point.roughness = point.roughness * point.roughness;
        point.roughness = clamp(point.roughness, 0.03f * 0.03f, 1.0f);
    }
    if (point.opacity > 0.999f) point.opacity = 1;
    if (voltransmission != zero3f) {
        point.voldensity = -log(clamp(voltransmission, 0.0001f, 1.0f)) /
                           volscale;
    } else {
        point.voldensity = zero3f;
    }
    return point;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

void merge_scene(yocto_scene& scene, const yocto_scene& merge) {
    auto offset_cameras      = scene.cameras.size();
    auto offset_textures     = scene.textures.size();
    auto offset_voltextures  = scene.voltextures.size();
    auto offset_materials    = scene.materials.size();
    auto offset_shapes       = scene.shapes.size();
    auto offset_subdivs      = scene.subdivs.size();
    auto offset_instances    = scene.instances.size();
    auto offset_environments = scene.environments.size();
    auto offset_nodes        = scene.nodes.size();
    auto offset_animations   = scene.animations.size();
    scene.cameras += merge.cameras;
    scene.textures += merge.textures;
    scene.voltextures += merge.voltextures;
    scene.materials += merge.materials;
    scene.shapes += merge.shapes;
    scene.subdivs += merge.subdivs;
    scene.instances += merge.instances;
    scene.environments += merge.environments;
    scene.nodes += merge.nodes;
    scene.animations += merge.animations;
    for (auto material_id = offset_materials;
         material_id < scene.materials.size(); material_id++) {
        auto& material = scene.materials[material_id];
        if (material.emission_texture >= 0)
            material.emission_texture += offset_textures;
        if (material.diffuse_texture >= 0)
            material.diffuse_texture += offset_textures;
        if (material.metallic_texture >= 0)
            material.metallic_texture += offset_textures;
        if (material.specular_texture >= 0)
            material.specular_texture += offset_textures;
        if (material.transmission_texture >= 0)
            material.transmission_texture += offset_textures;
        if (material.roughness_texture >= 0)
            material.roughness_texture += offset_textures;
        if (material.normal_texture >= 0)
            material.normal_texture += offset_textures;
        if (material.volume_density_texture >= 0)
            material.volume_density_texture += offset_voltextures;
    }
    for (auto subdiv_id = offset_subdivs; subdiv_id < scene.subdivs.size();
         subdiv_id++) {
        auto& subdiv = scene.subdivs[subdiv_id];
        if (subdiv.tesselated_shape >= 0)
            subdiv.tesselated_shape += offset_shapes;
        if (subdiv.displacement_texture >= 0)
            subdiv.displacement_texture += offset_textures;
    }
    for (auto instance_id = offset_instances;
         instance_id < scene.instances.size(); instance_id++) {
        auto& instance = scene.instances[instance_id];
        if (instance.shape >= 0) instance.shape += offset_shapes;
        if (instance.material >= 0) instance.material += offset_materials;
    }
    for (auto environment_id = offset_environments;
         environment_id < scene.environments.size(); environment_id++) {
        auto& environment = scene.environments[environment_id];
        if (environment.emission_texture >= 0)
            environment.emission_texture += offset_textures;
    }
    for (auto node_id = offset_nodes; node_id < scene.nodes.size(); node_id++) {
        auto& node = scene.nodes[node_id];
        if (node.parent >= 0) node.parent += offset_nodes;
        if (node.camera >= 0) node.camera += offset_cameras;
        if (node.instance >= 0) node.instance += offset_instances;
        if (node.environment >= 0) node.environment += offset_environments;
    }
    for (auto animation_id = offset_animations;
         animation_id < scene.animations.size(); animation_id++) {
        auto& animation = scene.animations[animation_id];
        for (auto& target : animation.node_targets)
            if (target >= 0) target += offset_nodes;
    }
}

string format_stats(const yocto_scene& scene, bool verbose) {
    auto stats = vector<pair<string, size_t>>{};
    stats += {"cameras", scene.cameras.size()};
    stats += {"shapes", scene.shapes.size()};
    stats += {"subdivs", scene.subdivs.size()};
    stats += {"instances", scene.instances.size()};
    stats += {"environments", scene.environments.size()};
    stats += {"textures", scene.textures.size()};
    stats += {"voltextures", scene.voltextures.size()};
    stats += {"materials", scene.materials.size()};
    stats += {"nodes", scene.nodes.size()};
    stats += {"animations", scene.animations.size()};

    auto accumulate = [](const auto& values, const auto& func) -> size_t {
        auto sum = (size_t)0;
        for (auto& value : values) sum += func(value);
        return sum;
    };

    stats += {"points", accumulate(scene.shapes,
                            [](auto& shape) { return shape.points.size(); })};
    stats += {"lines", accumulate(scene.shapes,
                           [](auto& shape) { return shape.lines.size(); })};
    stats += {"triangles", accumulate(scene.shapes, [](auto& shape) {
                  return shape.triangles.size();
              })};
    stats += {"quads", accumulate(scene.shapes,
                           [](auto& shape) { return shape.quads.size(); })};
    stats += {"fvquads", accumulate(scene.shapes, [](auto& shape) {
                  return shape.quads_positions.size();
              })};

    stats += {"texels4b", accumulate(scene.textures, [](auto& texture) {
                  return (size_t)texture.ldr_image.size().x *
                         (size_t)texture.ldr_image.size().x;
              })};
    stats += {"texels4f", accumulate(scene.textures, [](auto& texture) {
                  return (size_t)texture.hdr_image.size().x *
                         (size_t)texture.hdr_image.size().y;
              })};
    stats += {"voltexels", accumulate(scene.voltextures, [](auto& texture) {
                  return (size_t)texture.volume_data.size().x *
                         (size_t)texture.volume_data.size().y *
                         (size_t)texture.volume_data.size().z;
              })};

    auto str = ""s;
    for (auto& [key, value] : stats) {
        if (value == 0) continue;
        str += format("{:<15} {:>13n}\n", key + ":", value);
    }

    return str;
}

}  // namespace yocto
