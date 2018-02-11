//
// Implementation of Yocto/glTF. See header file for documentation.
//

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

#include "yocto_gltf.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR GLTF HIGH-LEVEL INTERFACE
// -----------------------------------------------------------------------------

namespace ygl {

// Math support
inline mat4f node_transform(const gltf_node* node) {
    return frame_to_mat(translation_frame(node->translation) * rotation_frame(node->rotation) *
           scaling_frame(node->scale)) * node->matrix;
}

// cleanup
gltf_mesh::~gltf_mesh() {
    for (auto e : shapes)
        if (e) delete e;
}

// cleanup
gltf_animation_group::~gltf_animation_group() {
    for (auto e : animations)
        if (e) delete e;
}

// cleanup
gltf_scene_group::~gltf_scene_group() {
    for (auto e : cameras)
        if (e) delete e;
    for (auto e : materials)
        if (e) delete e;
    for (auto e : meshes)
        if (e) delete e;
    for (auto e : textures)
        if (e) delete e;
    for (auto e : nodes)
        if (e) delete e;
    for (auto e : scenes)
        if (e) delete e;
    for (auto e : animations)
        if (e) delete e;
    for (auto e : skins)
        if (e) delete e;
}

// Flattens a gltf file into a flattened asset.
gltf_scene_group* gltf_to_scenes(const glTF* gltf, int scene_idx) {
    // clear asset
    auto scns = new gltf_scene_group();

    // convert images
    for (auto gtxt : gltf->images) {
        auto txt = new gltf_texture();
        txt->name = gtxt->name;
        txt->path = (startswith(gtxt->uri, "data:")) ? std::string("inlines") :
                                                       gtxt->uri;
        if (!gtxt->data.datab.empty()) {
            txt->ldr = image4b(gtxt->data.width, gtxt->data.height);
            for (auto j = 0; j < gtxt->data.height; j++) {
                for (auto i = 0; i < gtxt->data.width; i++) {
                    auto v = gtxt->data.datab.data() +
                             (gtxt->data.width * j + i) * gtxt->data.ncomp;
                    switch (gtxt->data.ncomp) {
                        case 1:
                            txt->ldr.at(i, j) = {v[0], v[0], v[0], 255};
                            break;
                        case 2: txt->ldr.at(i, j) = {v[0], v[1], 0, 255}; break;
                        case 3:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], 255};
                            break;
                        case 4:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        } else if (!gtxt->data.dataf.empty()) {
            txt->hdr = image4f(gtxt->data.width, gtxt->data.height);
            for (auto j = 0; j < gtxt->data.height; j++) {
                for (auto i = 0; i < gtxt->data.width; i++) {
                    auto v = gtxt->data.dataf.data() +
                             (gtxt->data.width * j + i) * gtxt->data.ncomp;
                    switch (gtxt->data.ncomp) {
                        case 1:
                            txt->hdr.at(i, j) = {v[0], v[0], v[0], 1};
                            break;
                        case 2: txt->hdr.at(i, j) = {v[0], v[1], 0, 1}; break;
                        case 3:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], 1};
                            break;
                        case 4:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        }
        scns->textures.push_back(txt);
    }

    // maps for translation
    static const auto filter_min_map =
        std::map<glTFSamplerMinFilter, gltf_texture_filter>{
            {glTFSamplerMinFilter::NotSet,
                gltf_texture_filter::linear_mipmap_linear},
            {glTFSamplerMinFilter::Linear, gltf_texture_filter::linear},
            {glTFSamplerMinFilter::Nearest, gltf_texture_filter::nearest},
            {glTFSamplerMinFilter::LinearMipmapLinear,
                gltf_texture_filter::linear_mipmap_linear},
            {glTFSamplerMinFilter::LinearMipmapNearest,
                gltf_texture_filter::linear_mipmap_nearest},
            {glTFSamplerMinFilter::NearestMipmapLinear,
                gltf_texture_filter::nearest_mipmap_linear},
            {glTFSamplerMinFilter::NearestMipmapNearest,
                gltf_texture_filter::nearest_mipmap_nearest},
        };
    static const auto filter_mag_map =
        std::map<glTFSamplerMagFilter, gltf_texture_filter>{
            {glTFSamplerMagFilter::NotSet, gltf_texture_filter::linear},
            {glTFSamplerMagFilter::Linear, gltf_texture_filter::linear},
            {glTFSamplerMagFilter::Nearest, gltf_texture_filter::nearest},
        };
    static const auto wrap_s_map =
        std::map<glTFSamplerWrapS, gltf_texture_wrap>{
            {glTFSamplerWrapS::NotSet, gltf_texture_wrap::repeat},
            {glTFSamplerWrapS::Repeat, gltf_texture_wrap::repeat},
            {glTFSamplerWrapS::ClampToEdge, gltf_texture_wrap::clamp},
            {glTFSamplerWrapS::MirroredRepeat, gltf_texture_wrap::mirror},
        };
    static const auto wrap_t_map =
        std::map<glTFSamplerWrapT, gltf_texture_wrap>{
            {glTFSamplerWrapT::NotSet, gltf_texture_wrap::repeat},
            {glTFSamplerWrapT::Repeat, gltf_texture_wrap::repeat},
            {glTFSamplerWrapT::ClampToEdge, gltf_texture_wrap::clamp},
            {glTFSamplerWrapT::MirroredRepeat, gltf_texture_wrap::mirror},
        };

    // add a texture
    auto add_texture = [gltf, scns](glTFTextureInfo* ginfo, gltf_texture*& txt,
                           gltf_texture_info*& info, bool normal = false,
                           bool occlusion = false) {
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt) return;
        txt = (!gtxt->source) ? nullptr : scns->textures[(int)gtxt->source];
        if (!txt) return;
        info = new gltf_texture_info();
        auto gsmp = gltf->get(gtxt->sampler);
        if (gsmp) {
            info->filter_mag = filter_mag_map.at(gsmp->magFilter);
            info->filter_min = filter_min_map.at(gsmp->minFilter);
            info->wrap_s = wrap_s_map.at(gsmp->wrapS);
            info->wrap_t = wrap_t_map.at(gsmp->wrapT);
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            info->scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            info->scale = ninfo->strength;
        }
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new gltf_material();
        mat->name = gmat->name;
        mat->emission = gmat->emissiveFactor;
        if (gmat->emissiveTexture) {
            add_texture(gmat->emissiveTexture, mat->emission_txt,
                mat->emission_txt_info);
        }
        if (gmat->pbrMetallicRoughness) {
            auto gmr = gmat->pbrMetallicRoughness;
            mat->metallic_roughness = new gltf_material_metallic_roughness();
            auto mr = mat->metallic_roughness;
            mr->base = {gmr->baseColorFactor[0], gmr->baseColorFactor[1],
                gmr->baseColorFactor[2]};
            mr->opacity = gmr->baseColorFactor[3];
            mr->metallic = gmr->metallicFactor;
            mr->roughness = gmr->roughnessFactor;
            if (gmr->baseColorTexture) {
                add_texture(
                    gmr->baseColorTexture, mr->base_txt, mr->base_txt_info);
            }
            if (gmr->metallicRoughnessTexture) {
                add_texture(gmr->metallicRoughnessTexture, mr->metallic_txt,
                    mr->metallic_txt_info);
            }
        }
        if (gmat->pbrSpecularGlossiness) {
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->specular_glossiness = new gltf_material_specular_glossiness();
            auto sg = mat->specular_glossiness;
            sg->diffuse = {gsg->diffuseFactor[0], gsg->diffuseFactor[1],
                gsg->diffuseFactor[2]};
            sg->opacity = gsg->diffuseFactor[3];
            sg->specular = gsg->specularFactor;
            sg->glossiness = gsg->glossinessFactor;
            if (gsg->diffuseTexture) {
                add_texture(
                    gsg->diffuseTexture, sg->diffuse_txt, sg->diffuse_txt_info);
            }
            if (gsg->specularGlossinessTexture) {
                add_texture(gsg->specularGlossinessTexture, sg->specular_txt,
                    sg->specular_txt_info);
            }
        }
        if (gmat->normalTexture) {
            add_texture(gmat->normalTexture, mat->normal_txt,
                mat->normal_txt_info, true, false);
        }
        if (gmat->occlusionTexture) {
            add_texture(gmat->occlusionTexture, mat->occlusion_txt,
                mat->occlusion_txt_info, false, true);
        }
        mat->double_sided = gmat->doubleSided;
        scns->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<gltf_shape*>>();
    for (auto gmesh : gltf->meshes) {
        auto msh = new gltf_mesh();
        // primitives
        for (auto gprim : gmesh->primitives) {
            auto prim = new gltf_shape();
            if (gprim->material) {
                prim->mat = scns->materials[(int)gprim->material];
            }
            // vertex data
            for (auto gattr : gprim->attributes) {
                auto semantic = gattr.first;
                auto vals = accessor_view(gltf, gltf->get(gattr.second));
                if (semantic == "POSITION") {
                    prim->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->pos.push_back(vals.getv3f(i));
                } else if (semantic == "NORMAL") {
                    prim->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->norm.push_back(vals.getv3f(i));
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    prim->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->texcoord.push_back(vals.getv2f(i));
                } else if (semantic == "TEXCOORD_1") {
                    prim->texcoord1.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->texcoord1.push_back(vals.getv2f(i));
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    prim->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->color.push_back(vals.getv4f(i, {0, 0, 0, 1}));
                } else if (semantic == "TANGENT") {
                    prim->tangsp.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->tangsp.push_back(vals.getv4f(i));
                } else if (semantic == "WEIGHTS_0") {
                    prim->skin_weights.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->skin_weights.push_back(vals.getv4f(i));
                } else if (semantic == "JOINTS_0") {
                    prim->skin_joints.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->skin_joints.push_back(vals.getv4i(i));
                } else if (semantic == "RADIUS") {
                    prim->radius.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->radius.push_back(vals.get(i, 0));
                } else {
                    // ignore
                }
            }
            // indices
            if (!gprim->indices) {
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        prim->triangles.reserve(prim->pos.size() / 3);
                        for (auto i = 0; i < prim->pos.size() / 3; i++) {
                            prim->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({0, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({i - 2, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        prim->lines.reserve(prim->pos.size() / 2);
                        for (auto i = 0; i < prim->pos.size() / 2; i++) {
                            prim->lines.push_back({i * 2 + 0, i * 2 + 1});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        prim->lines.reserve(prim->pos.size());
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                        prim->lines.back() = {(int)prim->pos.size() - 1, 0};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        prim->lines.reserve(prim->pos.size() - 1);
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        prim->points.reserve(prim->pos.size());
                        for (auto i = 0; i < prim->pos.size(); i++) {
                            prim->points.push_back(i);
                        }
                    } break;
                }
            } else {
                auto indices = accessor_view(gltf, gltf->get(gprim->indices));
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        prim->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++) {
                            prim->triangles.push_back({indices.geti(i * 3 + 0),
                                indices.geti(i * 3 + 1),
                                indices.geti(i * 3 + 2)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back({indices.geti(0),
                                indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back({indices.geti(i - 2),
                                indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        prim->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++) {
                            prim->lines.push_back({indices.geti(i * 2 + 0),
                                indices.geti(i * 2 + 1)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        prim->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back(
                                {indices.geti(i - 1), indices.geti(i)});
                        }
                        prim->lines.back() = {
                            indices.geti(indices.size() - 1), indices.geti(0)};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        prim->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back(
                                {indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        prim->points.reserve(indices.size());
                        for (auto i = 0; i < indices.size(); i++) {
                            prim->points.push_back(indices.geti(i));
                        }
                    } break;
                }
            }

            // morph targets
            int target_index = 0;
            for (auto& gtarget : gprim->targets) {
                auto target = new gltf_shape_morph();
                for (auto gattr : gtarget) {
                    auto semantic = gattr.first;
                    auto vals = accessor_view(gltf, gltf->get(gattr.second));
                    if (semantic == "POSITION") {
                        target->pos.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->pos.push_back(vals.getv3f(i));
                    } else if (semantic == "NORMAL") {
                        target->norm.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->norm.push_back(vals.getv3f(i));
                    } else if (semantic == "TANGENT") {
                        target->tangsp.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->tangsp.push_back(vals.getv3f(i));
                    } else {
                        // ignore
                    }
                }
                if (target_index < (int)gmesh->weights.size() - 1)
                    target->weight = gmesh->weights[target_index];
                target_index++;
                prim->morph_targets.push_back(target);
            }
            msh->shapes.push_back(prim);
        }
        scns->meshes.push_back(msh);
    }

    // convert cameras
    for (auto gcam : gltf->cameras) {
        auto cam = new gltf_camera();
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
        scns->cameras.push_back(cam);
    }

    // convert nodes
    for (auto gnode : gltf->nodes) {
        auto node = new gltf_node();
        node->name = gnode->name;
        node->cam =
            (!gnode->camera) ? nullptr : scns->cameras[(int)gnode->camera];
        node->msh = (!gnode->mesh) ? nullptr : scns->meshes[(int)gnode->mesh];
        node->translation = gnode->translation;
        node->rotation = gnode->rotation;
        node->scale = gnode->scale;
        node->matrix = gnode->matrix;
        node->morph_weights = gnode->weights;
        scns->nodes.push_back(node);
    }

    // set up children pointers
    auto is_root = std::vector<bool>(gltf->nodes.size(), true);
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnode = gltf->nodes[nid];
        auto node = scns->nodes[nid];
        for (auto n : gnode->children)
            node->children.push_back(scns->nodes[(int)n]);
        for (auto n : gnode->children) is_root[(int)n] = false;
    }

    // fix node morph weights
    for (auto node : scns->nodes) {
        if (!node->msh) continue;
        for (auto shp : node->msh->shapes) {
            if (node->morph_weights.size() < shp->morph_targets.size()) {
                node->morph_weights.resize(shp->morph_targets.size());
            }
        }
    }

    // convert animations
    for (auto ganim : gltf->animations) {
        auto anim_group = new gltf_animation_group();
        anim_group->name = ganim->name;
        std::map<vec2i, gltf_animation*> sampler_map;
        for (auto gchannel : ganim->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganim->get(gchannel->sampler);
                auto keyframes = new gltf_animation();
                auto input_view =
                    accessor_view(gltf, gltf->get(gsampler->input));
                keyframes->time.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    keyframes->time[i] = input_view.get(i);
                keyframes->interp =
                    (gltf_animation_interpolation)gsampler->interpolation;
                auto output_view =
                    accessor_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        keyframes->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->translation.push_back(
                                output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        keyframes->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->rotation.push_back(
                                (quat4f)output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        keyframes->scale.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->scale.push_back(output_view.getv3f(i));
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
                            keyframes->morph_weights.resize(
                                values.size() / ncomp);
                            for (auto i = 0;
                                 i < keyframes->morph_weights.size(); i++) {
                                keyframes->morph_weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    keyframes->morph_weights[i][j] =
                                        values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        // skip
                    }
                }
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = keyframes;
                anim_group->animations.push_back(keyframes);
            }
            sampler_map
                .at({(int)gchannel->sampler, (int)gchannel->target->path})
                ->nodes.push_back(scns->nodes[(int)gchannel->target->node]);
        }
        scns->animations.push_back(anim_group);
    }

    // convert skins
    for (auto gskin : gltf->skins) {
        auto skin = new gltf_skin();
        skin->name = gskin->name;
        for (auto gnode : gskin->joints)
            skin->joints.push_back(scns->nodes[(int)gnode]);
        skin->root = scns->nodes[(int)gskin->skeleton];
        if (!gskin->inverseBindMatrices) {
            skin->pose_matrices.assign(skin->joints.size(), identity_mat4f);
        } else {
            auto pose_matrix_view =
                accessor_view(gltf, gltf->get(gskin->inverseBindMatrices));
            skin->pose_matrices.resize(skin->joints.size());
            assert(pose_matrix_view.size() == skin->joints.size());
            assert(pose_matrix_view.ncomp() == 16);
            for (auto i = 0; i < pose_matrix_view.size(); i++) {
                skin->pose_matrices[i] = pose_matrix_view.getm4f(i);
            }
        }
        scns->skins.push_back(skin);
    }

    // set skin pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        if (!gltf->nodes[nid]->skin) continue;
        scns->nodes[nid]->skn = scns->skins[(int)gltf->nodes[nid]->skin];
    }

    // convert scenes
    for (auto gscn : gltf->scenes) {
        auto scn = new gltf_scene();
        scn->name = gscn->name;
        for (auto n : gscn->nodes) scn->nodes.push_back(scns->nodes[(int)n]);
        scns->scenes.push_back(scn);
    }
    if (gltf->scene) { scns->default_scene = scns->scenes[(int)gltf->scene]; }

    // update transforms
    update_transforms(scns);

    return scns;
}

// helper
template <typename T>
static inline int index(const std::vector<T*>& vec, T* val) {
    auto pos = std::find(vec.begin(), vec.end(), val);
    if (pos == vec.end()) return -1;
    return (int)(pos - vec.begin());
}

// Unflattnes gltf
glTF* scenes_to_gltf(const gltf_scene_group* scns,
    const std::string& buffer_uri, bool separate_buffers) {
    auto gltf = std::unique_ptr<glTF>(new glTF());

    // add asset info
    gltf->asset = new glTFAsset();
    gltf->asset->generator = "Yocto/gltf";
    gltf->asset->version = "2.0";

    // convert cameras
    for (auto cam : scns->cameras) {
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
    for (auto txt : scns->textures) {
        auto gimg = new glTFImage();
        gimg->uri = txt->path;
        if (txt->hdr) {
            gimg->data.width = txt->hdr.width();
            gimg->data.height = txt->hdr.height();
            gimg->data.ncomp = 4;
            gimg->data.dataf.assign((uint8_t*)data(txt->hdr),
                (uint8_t*)data(txt->hdr) +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (txt->ldr) {
            gimg->data.width = txt->ldr.width();
            gimg->data.height = txt->ldr.height();
            gimg->data.ncomp = 4;
            gimg->data.datab.assign((uint8_t*)data(txt->ldr),
                (uint8_t*)data(txt->ldr) +
                    txt->ldr.width() * txt->ldr.height() * 4);
        }
        gltf->images.push_back(gimg);
    }

    // conversion maps
    static const auto wrap_s_map =
        std::map<gltf_texture_wrap, glTFSamplerWrapS>{
            {gltf_texture_wrap::repeat, glTFSamplerWrapS::Repeat},
            {gltf_texture_wrap::clamp, glTFSamplerWrapS::ClampToEdge},
            {gltf_texture_wrap::mirror, glTFSamplerWrapS::MirroredRepeat},
        };
    static const auto wrap_t_map =
        std::map<gltf_texture_wrap, glTFSamplerWrapT>{
            {gltf_texture_wrap::repeat, glTFSamplerWrapT::Repeat},
            {gltf_texture_wrap::clamp, glTFSamplerWrapT::ClampToEdge},
            {gltf_texture_wrap::mirror, glTFSamplerWrapT::MirroredRepeat},
        };
    static const auto texture_min_map =
        std::map<gltf_texture_filter, glTFSamplerMinFilter>{
            {gltf_texture_filter::linear, glTFSamplerMinFilter::Linear},
            {gltf_texture_filter::nearest, glTFSamplerMinFilter::Nearest},
            {gltf_texture_filter::linear_mipmap_linear,
                glTFSamplerMinFilter::LinearMipmapLinear},
            {gltf_texture_filter::linear_mipmap_nearest,
                glTFSamplerMinFilter::LinearMipmapNearest},
            {gltf_texture_filter::nearest_mipmap_linear,
                glTFSamplerMinFilter::NearestMipmapNearest},
            {gltf_texture_filter::nearest_mipmap_nearest,
                glTFSamplerMinFilter::NearestMipmapNearest},
        };
    static const auto texture_mag_map =
        std::map<gltf_texture_filter, glTFSamplerMagFilter>{
            {gltf_texture_filter::linear, glTFSamplerMagFilter::Linear},
            {gltf_texture_filter::nearest, glTFSamplerMagFilter::Nearest},
        };

    // add a texture and sampler
    auto add_texture = [&gltf, scns](
                           gltf_texture* txt, gltf_texture_info* info) {
        if (!txt) return glTFid<glTFTexture>{-1};
        auto gtxt = new glTFTexture();
        gtxt->name = txt->name;
        gtxt->source = glTFid<glTFImage>(index(scns->textures, txt));
        if (info && !info->is_default()) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = wrap_s_map.at(info->wrap_s);
            gsmp->wrapT = wrap_t_map.at(info->wrap_t);
            gsmp->minFilter = texture_min_map.at(info->filter_min);
            gsmp->magFilter = texture_mag_map.at(info->filter_mag);
            gtxt->sampler = glTFid<glTFSampler>((int)gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        return glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
    };

    // convert materials
    for (auto mat : scns->materials) {
        auto gmat = new glTFMaterial();
        gmat->name = mat->name;
        gmat->emissiveFactor = mat->emission;
        if (mat->emission_txt) {
            gmat->emissiveTexture = new glTFTextureInfo();
            gmat->emissiveTexture->index =
                add_texture(mat->emission_txt, mat->emission_txt_info);
        }
        if (mat->metallic_roughness) {
            gmat->pbrMetallicRoughness = new glTFMaterialPbrMetallicRoughness();
            auto gmr = gmat->pbrMetallicRoughness;
            auto mr = mat->metallic_roughness;
            gmr->baseColorFactor = {
                mr->base[0], mr->base[1], mr->base[2], mr->opacity};
            gmr->metallicFactor = mr->metallic;
            gmr->roughnessFactor = mr->roughness;
            if (mr->base_txt) {
                gmr->baseColorTexture = new glTFTextureInfo();
                gmr->baseColorTexture->index =
                    add_texture(mr->base_txt, mr->base_txt_info);
            }
            if (mr->metallic_txt) {
                gmr->metallicRoughnessTexture = new glTFTextureInfo();
                gmr->metallicRoughnessTexture->index =
                    add_texture(mr->metallic_txt, mr->metallic_txt_info);
            }
        }
        if (mat->specular_glossiness) {
            gmat->pbrSpecularGlossiness =
                new glTFMaterialPbrSpecularGlossiness();
            auto gsg = gmat->pbrSpecularGlossiness;
            auto sg = mat->specular_glossiness;
            gsg->diffuseFactor = {
                sg->diffuse[0], sg->diffuse[1], sg->diffuse[2], sg->opacity};
            gsg->specularFactor = sg->specular;
            gsg->glossinessFactor = sg->glossiness;
            if (sg->diffuse_txt) {
                gsg->diffuseTexture = new glTFTextureInfo();
                gsg->diffuseTexture->index =
                    add_texture(sg->diffuse_txt, sg->diffuse_txt_info);
            }
            if (sg->specular_txt) {
                gsg->specularGlossinessTexture = new glTFTextureInfo();
                gsg->specularGlossinessTexture->index =
                    add_texture(sg->specular_txt, sg->specular_txt_info);
            }
        }
        if (mat->normal_txt) {
            gmat->normalTexture = new glTFMaterialNormalTextureInfo();
            gmat->normalTexture->index =
                add_texture(mat->normal_txt, mat->normal_txt_info);
            if (gmat->normalTexture && mat->normal_txt_info)
                gmat->normalTexture->scale = mat->normal_txt_info->scale;
        }
        if (mat->occlusion_txt) {
            gmat->occlusionTexture = new glTFMaterialOcclusionTextureInfo();
            gmat->occlusionTexture->index =
                add_texture(mat->occlusion_txt, mat->occlusion_txt_info);
            if (gmat->occlusionTexture && mat->occlusion_txt_info)
                gmat->occlusionTexture->strength =
                    mat->occlusion_txt_info->scale;
        }
        gmat->doubleSided = mat->double_sided;
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
    auto add_accessor = [&gltf](glTFBuffer* gbuffer, const std::string& name,
                            glTFAccessorType type,
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
            switch (type) {
                case glTFAccessorType::Scalar: {
                    auto bbox = make_bbox(count, (vec1f*)data);
                    accessor->min = {bbox.min.x};
                    accessor->max = {bbox.max.x};
                } break;
                case glTFAccessorType::Vec2: {
                    auto bbox = make_bbox(count, (vec2f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y};
                    accessor->max = {bbox.max.x, bbox.max.y};
                } break;
                case glTFAccessorType::Vec3: {
                    auto bbox = make_bbox(count, (vec3f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y, bbox.min.z};
                    accessor->max = {bbox.max.x, bbox.max.y, bbox.max.z};
                } break;
                case glTFAccessorType::Vec4: {
                    auto bbox = make_bbox(count, (vec4f*)data);
                    accessor->min = {
                        bbox.min.x, bbox.min.y, bbox.min.z, bbox.min.w};
                    accessor->max = {
                        bbox.max.x, bbox.max.y, bbox.max.z, bbox.max.w};
                } break;
                default: break;
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // convert meshes
    for (auto msh : scns->meshes) {
        auto gbuffer = add_opt_buffer(msh->path);
        auto gmesh = new glTFMesh();
        gmesh->name = msh->name;
        for (auto j = 0; j < msh->shapes.size(); j++) {
            auto gprim = msh->shapes[j];
            auto pid = msh->name + "_" + std::to_string(j);
            auto prim = new glTFMeshPrimitive();
            prim->material =
                glTFid<glTFMaterial>(index(scns->materials, gprim->mat));
            if (!gprim->pos.empty())
                prim->attributes["POSITION"] = add_accessor(gbuffer,
                    pid + "_pos", glTFAccessorType::Vec3,
                    glTFAccessorComponentType::Float, (int)gprim->pos.size(),
                    sizeof(vec3f), gprim->pos.data(), true);
            if (!gprim->norm.empty())
                prim->attributes["NORMAL"] = add_accessor(gbuffer,
                    pid + "_norm", glTFAccessorType::Vec3,
                    glTFAccessorComponentType::Float, (int)gprim->norm.size(),
                    sizeof(vec3f), gprim->norm.data(), false);
            if (!gprim->texcoord.empty())
                prim->attributes["TEXCOORD_0"] = add_accessor(gbuffer,
                    pid + "_texcoord", glTFAccessorType::Vec2,
                    glTFAccessorComponentType::Float,
                    (int)gprim->texcoord.size(), sizeof(vec2f),
                    gprim->texcoord.data(), false);
            if (!gprim->texcoord1.empty())
                prim->attributes["TEXCOORD_1"] = add_accessor(gbuffer,
                    pid + "_texcoord1", glTFAccessorType::Vec2,
                    glTFAccessorComponentType::Float,
                    (int)gprim->texcoord1.size(), sizeof(vec2f),
                    gprim->texcoord1.data(), false);
            if (!gprim->color.empty())
                prim->attributes["COLOR_0"] = add_accessor(gbuffer,
                    pid + "_color", glTFAccessorType::Vec4,
                    glTFAccessorComponentType::Float, (int)gprim->color.size(),
                    sizeof(vec4f), gprim->color.data(), false);
            if (!gprim->skin_weights.empty())
                prim->attributes["WEIGHTS_0"] = add_accessor(gbuffer,
                    pid + "_skin_weights", glTFAccessorType::Vec4,
                    glTFAccessorComponentType::Float,
                    (int)gprim->skin_weights.size(), sizeof(vec4f),
                    gprim->skin_weights.data(), false);
            if (!gprim->skin_joints.empty()) {
                using ushort = unsigned short;
                auto joints_short = std::vector<std::array<ushort, 4>>();
                joints_short.reserve(gprim->skin_joints.size());
                for (auto&& j : gprim->skin_joints)
                    joints_short.push_back(
                        {(ushort)j.x, (ushort)j.y, (ushort)j.z, (ushort)j.w});
                prim->attributes["JOINTS_0"] = add_accessor(gbuffer,
                    pid + "_skin_joints", glTFAccessorType::Vec4,
                    glTFAccessorComponentType::UnsignedShort,
                    (int)joints_short.size(), sizeof(ushort) * 4,
                    joints_short.data(), false);
            }
            if (!gprim->radius.empty())
                prim->attributes["RADIUS"] = add_accessor(gbuffer,
                    pid + "_radius", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)gprim->radius.size(),
                    sizeof(float), gprim->radius.data(), false);
            // auto elem_as_uint = gprim->pos.size() >
            // std::numeric_limits<unsigned short>::max();
            if (!gprim->points.empty()) {
                prim->indices = add_accessor(gbuffer, pid + "_points",
                    glTFAccessorType::Scalar,
                    glTFAccessorComponentType::UnsignedInt,
                    (int)gprim->points.size(), sizeof(int),
                    (int*)gprim->points.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Points;
            } else if (!gprim->lines.empty()) {
                prim->indices = add_accessor(gbuffer, pid + "_lines",
                    glTFAccessorType::Scalar,
                    glTFAccessorComponentType::UnsignedInt,
                    (int)gprim->lines.size() * 2, sizeof(int),
                    (int*)gprim->lines.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Lines;
            } else if (!gprim->triangles.empty()) {
                prim->indices = add_accessor(gbuffer, pid + "_triangles",
                    glTFAccessorType::Scalar,
                    glTFAccessorComponentType::UnsignedInt,
                    (int)gprim->triangles.size() * 3, sizeof(int),
                    (int*)gprim->triangles.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Triangles;
            } else {
                assert(false);
            }
            auto target_index = 0;
            for (auto target : gprim->morph_targets) {
                auto tid = std::to_string(target_index++);
                prim->targets.push_back({});
                if (!target->pos.empty()) {
                    prim->targets.back()["POSITION"] = add_accessor(gbuffer,
                        pid + "_" + tid + "_pos", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->pos.size(), sizeof(vec3f),
                        target->pos.data(), true);
                } else if (!target->norm.empty()) {
                    prim->targets.back()["NORMAL"] = add_accessor(gbuffer,
                        pid + "_" + tid + "_norm", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->norm.size(), sizeof(vec3f),
                        target->norm.data(), true);
                } else if (!target->tangsp.empty()) {
                    prim->targets.back()["TANGENT"] = add_accessor(gbuffer,
                        pid + "_" + tid + "_tang", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->tangsp.size(), sizeof(vec3f),
                        target->tangsp.data(), true);
                } else {
                    assert(false);
                }
            }
            gmesh->primitives.push_back(prim);
        }
        gltf->meshes.push_back(gmesh);
    }

    // nodes
    for (auto node : scns->nodes) {
        auto gnode = new glTFNode();
        gnode->name = node->name;
        gnode->camera = glTFid<glTFCamera>(index(scns->cameras, node->cam));
        gnode->mesh = glTFid<glTFMesh>(index(scns->meshes, node->msh));
        gnode->matrix = node->matrix;
        gnode->translation = node->translation;
        gnode->rotation = node->rotation;
        gnode->scale = node->scale;
        gnode->weights = node->morph_weights;
        gltf->nodes.push_back(gnode);
    }

    // fix children
    auto nid = 0;
    for (auto node : scns->nodes) {
        auto gnode = gltf->nodes[nid++];
        for (auto child : node->children)
            gnode->children.push_back(
                glTFid<glTFNode>(index(scns->nodes, child)));
    }

    // skins
    for (auto sk : scns->skins) {
        auto gsk = new glTFSkin();
        auto gbuffer = add_opt_buffer(sk->path);
        gsk->name = sk->name;
        gsk->skeleton = glTFid<glTFNode>(index(scns->nodes, sk->root));
        for (auto joint : sk->joints) {
            gsk->joints.push_back(glTFid<glTFNode>(index(scns->nodes, joint)));
        }
        if (!sk->pose_matrices.empty()) {
            gsk->inverseBindMatrices = add_accessor(gbuffer, sk->name + "_pose",
                glTFAccessorType::Mat4, glTFAccessorComponentType::Float,
                (int)sk->pose_matrices.size(), sizeof(mat4f),
                sk->pose_matrices.data(), false);
        }
        gltf->skins.push_back(gsk);
    }

    // fix skin references
    nid = 0;
    for (auto node : scns->nodes) {
        auto gnode = gltf->nodes[nid++];
        if (node->skn) {
            gnode->skin = glTFid<glTFSkin>{index(scns->skins, node->skn)};
        }
    }

    // interpolation map
    static const auto interpolation_map = std::map<gltf_animation_interpolation,
        glTFAnimationSamplerInterpolation>{
        {gltf_animation_interpolation::step,
            glTFAnimationSamplerInterpolation::Step},
        {gltf_animation_interpolation::linear,
            glTFAnimationSamplerInterpolation::Linear},
        {gltf_animation_interpolation::cubic,
            glTFAnimationSamplerInterpolation::CubicSpline},
        {gltf_animation_interpolation::catmull_rom,
            glTFAnimationSamplerInterpolation::CatmullRomSpline},
    };

    // animation
    for (auto anims : scns->animations) {
        auto ganim = new glTFAnimation();
        ganim->name = anims->name;
        auto gbuffer = add_opt_buffer(anims->path);
        auto count = 0;
        for (auto anim : anims->animations) {
            auto aid = ganim->name + "_" + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input =
                add_accessor(gbuffer, aid + "_time", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)anim->time.size(),
                    sizeof(float), anim->time.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!anim->translation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anim->translation.size(), sizeof(vec3f),
                    anim->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!anim->rotation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)anim->rotation.size(), sizeof(vec4f),
                    anim->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!anim->scale.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anim->scale.size(), sizeof(vec3f), anim->scale.data(),
                    false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!anim->morph_weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(
                    anim->morph_weights.size() * anim->morph_weights[0].size());
                for (auto i = 0; i < anim->morph_weights.size(); i++) {
                    values.insert(values.end(), anim->morph_weights[i].begin(),
                        anim->morph_weights[i].end());
                }
                gsmp->output = add_accessor(gbuffer, aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
            }
            gsmp->interpolation = interpolation_map.at(anim->interp);
            for (auto node : anim->nodes) {
                auto gchan = new glTFAnimationChannel();
                gchan->sampler =
                    glTFid<glTFAnimationSampler>{(int)ganim->samplers.size()};
                gchan->target = new glTFAnimationChannelTarget();
                gchan->target->node =
                    glTFid<glTFNode>{index(scns->nodes, node)};
                gchan->target->path = path;
                ganim->channels.push_back(gchan);
            }
            ganim->samplers.push_back(gsmp);
        }
        gltf->animations.push_back(ganim);
    }

    // scenes
    for (auto scn : scns->scenes) {
        auto gscn = new glTFScene();
        gscn->name = scn->name;
        for (auto child : scn->nodes)
            gscn->nodes.push_back(glTFid<glTFNode>(index(scns->nodes, child)));
        gltf->scenes.push_back(gscn);
    }
    gltf->scene = glTFid<glTFScene>(index(scns->scenes, scns->default_scene));

    // done
    return gltf.release();
}

// Load scene
gltf_scene_group* load_scenes(
    const std::string& filename, bool load_textures, bool skip_missing) {
    auto ext = path_extension(filename);
    auto gltf = std::unique_ptr<glTF>();
    if (ext != ".glb") {
        gltf = std::unique_ptr<glTF>(
            load_gltf(filename, true, load_textures, skip_missing));
    } else {
        gltf = std::unique_ptr<glTF>(
            load_binary_gltf(filename, true, load_textures, skip_missing));
    }
    if (!gltf) return nullptr;
    return gltf_to_scenes(gltf.get());
}

/// Save scene
void save_scenes(const std::string& filename, const std::string& buffer_uri,
    const gltf_scene_group* scn, bool save_textures, bool separate_buffers) {
    auto gltf = std::unique_ptr<glTF>(
        scenes_to_gltf(scn, buffer_uri, separate_buffers));
    save_gltf(filename, gltf.get(), true, save_textures);
}

// Computes a scene bounding box
bbox3f compute_scene_bounds(const gltf_scene_group* scn) {
    auto bbox_meshes = std::map<gltf_mesh*, bbox3f>{{nullptr, invalid_bbox3f}};
    for (auto mesh : scn->meshes) {
        bbox_meshes[mesh] = invalid_bbox3f;
        auto& bbox = bbox_meshes[mesh];
        for (auto shp : mesh->shapes)
            for (auto& p : shp->pos) bbox += p;
    }
    auto bbox = invalid_bbox3f;
    if (!scn->nodes.empty()) {
        for (auto ist : scn->nodes) {
            if (ist->msh)
                bbox += transform_bbox(ist->xform(), bbox_meshes.at(ist->msh));
        }
    } else {
        for (auto mesh : scn->meshes) bbox += bbox_meshes[mesh];
    }
    return bbox;
}

// Add missing data to the scene.
void add_normals(gltf_scene_group* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (!shp->norm.empty()) continue;
            shp->norm.resize(shp->pos.size(), {0, 0, 1});
            if (!shp->lines.empty() || !shp->triangles.empty()) {
                shp->norm =
                    compute_normals(shp->lines, shp->triangles, {}, shp->pos);
            }
        }
    }
}

// Add missing data to the scene.
void add_tangent_space(gltf_scene_group* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (!shp->mat) continue;
            if (shp->triangles.empty()) continue;
            if (!shp->tangsp.empty() || shp->texcoord.empty() ||
                !shp->mat->normal_txt)
                continue;
            shp->tangsp.resize(shp->pos.size());
            shp->tangsp = compute_tangent_frames(
                shp->triangles, shp->pos, shp->norm, shp->texcoord);
        }
    }
}

// Add missing data to the scene.
void add_radius(gltf_scene_group* scn, float radius) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (shp->points.empty() && shp->lines.empty()) continue;
            if (!shp->radius.empty()) continue;
            shp->radius.resize(shp->pos.size(), radius);
        }
    }
}

// Add missing data to the scene.
void add_texture_data(gltf_scene_group* scn) {
    for (auto txt : scn->textures) {
        if (!txt->hdr && !txt->ldr) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->ldr = image4b(1, 1, {255, 255, 255, 255});
        }
    }
}

// Add missing data to the scene.
void add_nodes(gltf_scene_group* scn) {
    if (!scn->nodes.empty()) return;
    for (auto mesh : scn->meshes) {
        auto ist = new gltf_node();
        ist->name = mesh->name;
        ist->msh = mesh;
        scn->nodes.push_back(ist);
    }
}

// Add missing data to the scene.
void add_scene(gltf_scene_group* scn) {
    if (!scn->scenes.empty()) return;
    auto s = new gltf_scene();
    s->name = "scene";
    update_transforms(scn);
    for (auto n : scn->nodes) {
        if (!n->parent) s->nodes.push_back(n);
    }
}

// Add missing data to the scene.
void add_names(gltf_scene_group* scn) {
    auto cid = 0;
    for (auto cam : scn->cameras) {
        if (cam->name.empty())
            cam->name = "<camera " + std::to_string(cid) + ">";
        cid++;
    }

    auto tid = 0;
    for (auto texture : scn->textures) {
        if (texture->name.empty())
            texture->name = "<texture " + std::to_string(tid) + ">";
        tid++;
    }

    auto mid = 0;
    for (auto mat : scn->materials) {
        if (mat->name.empty())
            mat->name = "<material " + std::to_string(mid) + ">";
        mid++;
    }

    auto mmid = 0;
    for (auto mesh : scn->meshes) {
        if (mesh->name.empty())
            mesh->name = "<mesh " + std::to_string(mmid) + ">";
        mmid++;
        auto sid = 0;
        for (auto shp : mesh->shapes) {
            if (shp->name.empty())
                shp->name = "<shape " + std::to_string(sid) + ">";
            sid++;
        }
    }

    auto nid = 0;
    for (auto node : scn->nodes) {
        if (node->name.empty())
            node->name = "<node " + std::to_string(nid) + ">";
        nid++;
    }

    auto sid = 0;
    for (auto scene : scn->scenes) {
        if (scene->name.empty())
            scene->name = "<scene " + std::to_string(sid) + ">";
        sid++;
    }
}

// Add a default camera that views the entire scene.
void add_default_cameras(gltf_scene_group* scns) {
    for (auto scn : scns->scenes) {
        auto cams = get_camera_nodes(scn);
        if (cams.empty()) {
            // TODO: scene bounds
            auto bbox = bbox3f{compute_scene_bounds(scns)};
            auto bbox_center = (bbox.max + bbox.min) / 2.0f;
            auto bbox_size = bbox.max - bbox.min;
            auto bbox_msize =
                max(bbox_size[0], max(bbox_size[1], bbox_size[2]));
            // set up camera
            auto cam = new gltf_camera();
            cam->name = scn->name + " camera";
            auto camera_dir = vec3f{1, 0.4f, 1};
            auto from = camera_dir * bbox_msize + bbox_center;
            auto to = bbox_center;
            auto up = vec3f{0, 1, 0};
            cam->ortho = false;
            cam->aspect = 16.0f / 9.0f;
            cam->yfov = 2 * atanf(0.5f);
            cam->aperture = 0;
            cam->focus = length(to - from);
            auto node = new gltf_node();
            node->matrix = frame_to_mat(lookat_frame(from, to, up));
            node->cam = cam;
            node->name = cam->name;
            scns->cameras.push_back(cam);
            scns->nodes.push_back(node);
            scn->nodes.push_back(node);
        }
    }
}

// Update node hierarchy
void update_node_hierarchy(gltf_scene_group* scns) {
    for (auto node : scns->nodes) node->parent = nullptr;
    for (auto node : scns->nodes) {
        for (auto child : node->children) child->parent = node;
    }
}

// Update node trasforms
void update_transforms(gltf_node* ist) {
    ist->_local_xform = node_transform(ist);
    if (ist->parent) {
        ist->_xform = ist->parent->_xform * ist->_local_xform;
    } else {
        ist->_xform = ist->_local_xform;
    }
    for (auto child : ist->children) update_transforms(child);
}

// Update node trasforms
void update_transforms(gltf_scene_group* scns) {
    update_node_hierarchy(scns);
    for (auto node : scns->nodes) { update_transforms(node); }
}

// Evalute interpolated values
void update_animated_node_transforms(const gltf_animation* anim, float time) {
    time = clamp(time, anim->time.front(), anim->time.back() - 0.001f);
    // get time slice
    auto i1 = 0, i2 = 0;
    auto t = 0.0f;
    auto interp = anim->interp;
    if (time <= anim->time.front()) {
        interp = gltf_animation_interpolation::step;
    } else if (time >= anim->time.back()) {
        i1 = (int)anim->time.size() - 1;
        i2 = (int)anim->time.size() - 2;
        interp = gltf_animation_interpolation::step;
    } else {
        for (i2 = 0; i2 < anim->time.size() && anim->time[i2] < time; i2++)
            ;
        i1 = i2 - 1;
        t = (time - anim->time[i1]) / (anim->time[i2] - anim->time[i1]);
    }

    // apply transforms
    if (!anim->translation.empty()) {
        auto trans = vec3f{0, 0, 0};
        switch (interp) {
            case gltf_animation_interpolation::step: {
                trans = anim->translation[i1];
            } break;
            case gltf_animation_interpolation::linear: {
                trans =
                    anim->translation[i1] * (1 - t) + anim->translation[i2] * t;
            } break;
            case gltf_animation_interpolation::catmull_rom: {
            } break;
            case gltf_animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->translation = trans;
    } else if (!anim->rotation.empty()) {
        auto rot = quat4f{0, 0, 0, 1};
        switch (interp) {
            case gltf_animation_interpolation::step: {
                rot = anim->rotation[i1];
            } break;
            case gltf_animation_interpolation::linear: {
                rot = slerp(anim->rotation[i1], anim->rotation[i2], t);
            } break;
            case gltf_animation_interpolation::catmull_rom: {
            } break;
            case gltf_animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->rotation = rot;
    } else if (!anim->scale.empty()) {
        auto scale = vec3f{1, 1, 1};
        switch (interp) {
            case gltf_animation_interpolation::step: {
                scale = anim->scale[i1];
            } break;
            case gltf_animation_interpolation::linear: {
                scale = anim->scale[i1] * (1 - t) + anim->scale[i2] * t;
            } break;
            case gltf_animation_interpolation::catmull_rom: {
            } break;
            case gltf_animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->scale = scale;
    } else if (!anim->morph_weights.empty()) {
        auto weights = std::vector<float>(anim->morph_weights[0].size());
        switch (interp) {
            case gltf_animation_interpolation::step: {
                weights = anim->morph_weights[i1];
            } break;
            case gltf_animation_interpolation::linear: {
                for (auto i = 0; i < weights.size(); i++) {
                    weights[i] = anim->morph_weights[i1][i] * (1 - t) +
                                 anim->morph_weights[i2][i] * t;
                }
            } break;
            case gltf_animation_interpolation::catmull_rom: {
            } break;
            case gltf_animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->morph_weights = weights;
    } else {
    }
}

// Compute animation
void update_animated_transforms(gltf_scene_group* scns, float time) {
    for (auto anim_group : scns->animations) {
        for (auto anim : anim_group->animations) {
            update_animated_node_transforms(anim, time);
        }
    }
}

// Update skin trasforms
void update_skin_transforms(gltf_node* ist, gltf_node* parent) {
    ist->_skin_xform = node_transform(ist);
    if (parent) ist->_skin_xform = parent->_skin_xform * ist->_skin_xform;
    for (auto child : ist->children) update_skin_transforms(child, ist);
}

// Skin transforms (local-to-object)
std::vector<mat4f> get_skin_transforms(
    const gltf_skin* sk, const mat4f& xform) {
    auto ret = std::vector<mat4f>(sk->joints.size());
    update_skin_transforms(sk->root, nullptr);
    auto inv_root = inverse(xform);
    for (auto i = 0; i < sk->joints.size(); i++) {
        if (!sk->pose_matrices.empty()) {
            ret[i] = inv_root * sk->joints[i]->_xform * sk->pose_matrices[i];
        } else {
            ret[i] = inv_root * sk->joints[i]->_xform;
        }
    }
    return ret;
}

// Compute shape morphing
void compute_morphing_deformation(const gltf_shape* shp,
    const std::vector<float>& weights, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec4f>& tangsp) {
    pos = shp->pos;
    norm = shp->norm;
    tangsp = shp->tangsp;
    for (auto idx = 0; idx < shp->morph_targets.size(); idx++) {
        auto morph = shp->morph_targets[idx];
        auto weight = (idx < weights.size()) ? weights[idx] : morph->weight;
        if (weight == 0) continue;
        if (!morph->pos.empty()) {
            for (auto i = 0; i < pos.size(); i++) {
                pos[i] += weight * morph->pos[i];
            }
        }
        if (!morph->norm.empty()) {
            for (auto i = 0; i < pos.size(); i++) {
                norm[i] += weight * morph->norm[i];
            }
        }
        if (!morph->tangsp.empty()) {
            for (auto i = 0; i < tangsp.size(); i++) {
                *(vec3f*)(&tangsp[i]) += weight * morph->tangsp[i];
            }
        }
    }
}

// Animation times
vec2f get_animation_bounds(const gltf_scene_group* scns) {
    auto range = vec2f{0, 0};
    for (auto anim_group : scns->animations) {
        for (auto anim : anim_group->animations) {
            range[0] = std::min(anim->time.front(), range[0]);
            range[1] = std::max(anim->time.back(), range[1]);
        }
    }
    return range;
}

// Helpder to filter nodes
void get_nodes(const std::vector<gltf_node*> nodes, bool has_camera,
    bool has_mesh, std::vector<gltf_node*>* filtered) {
    for (auto node : nodes) {
        if (has_camera && node->cam) filtered->push_back(node);
        if (has_mesh && node->msh) filtered->push_back(node);
        get_nodes(node->children, has_camera, has_mesh, filtered);
    }
}

// Helpder to filter nodes
std::vector<gltf_node*> get_nodes(
    const std::vector<gltf_node*> nodes, bool has_camera, bool has_mesh) {
    std::vector<gltf_node*> filtered;
    get_nodes(nodes, has_camera, has_mesh, &filtered);
    return filtered;
}

// Get a list of nodes with meshes
std::vector<gltf_node*> get_mesh_nodes(const gltf_scene* scn) {
    if (!scn) return {};
    return get_nodes(scn->nodes, false, true);
}

// Get a list of nodes with cameras
std::vector<gltf_node*> get_camera_nodes(const gltf_scene* scn) {
    if (!scn) return {};
    return get_nodes(scn->nodes, true, false);
}

// Set unique path names for outputting separate buffers
void add_unique_path_names(
    gltf_scene_group* scns, const std::string& buffer_uri) {
    auto mid = 0;
    for (auto msh : scns->meshes) {
        msh->path = buffer_uri + "mesh_" + std::to_string(mid++) + "_" +
                    msh->name + ".bin";
    }
    auto aid = 0;
    for (auto anm : scns->animations) {
        anm->path = buffer_uri + "anim_" + std::to_string(aid++) + "_" +
                    anm->name + ".bin";
    }
    auto sid = 0;
    for (auto skn : scns->meshes) {
        skn->path = buffer_uri + "skin_" + std::to_string(sid++) + "_" +
                    skn->name + ".bin";
    }
}

// Convert materials to spec gloss
void add_spec_gloss(gltf_scene_group* scns) {
    auto txts = std::map<std::pair<gltf_texture*, gltf_texture*>,
        std::pair<gltf_texture*, gltf_texture*>>();
    for (auto mat : scns->materials) {
        if (mat->specular_glossiness) continue;
        if (!mat->metallic_roughness) continue;
        mat->specular_glossiness = new gltf_material_specular_glossiness();
        auto mr = mat->metallic_roughness;
        auto sg = mat->specular_glossiness;
        if (mr->base_txt || mr->metallic_txt) {
            sg->diffuse = {1, 1, 1};
            sg->opacity = 1;
            sg->specular = {1, 1, 1};
            sg->glossiness = 1;
            auto mr_txt = std::pair<gltf_texture*, gltf_texture*>{
                mr->base_txt, mr->metallic_txt};
            if (txts.find(mr_txt) == txts.end()) {
                auto w = 0, h = 0;
                if (mr->base_txt) {
                    w = max(w, mr->base_txt->width());
                    h = max(h, mr->base_txt->height());
                }
                if (mr->metallic_txt) {
                    w = max(w, mr->metallic_txt->width());
                    h = max(h, mr->metallic_txt->height());
                }
                auto diff = new gltf_texture();
                diff->ldr = image4b(w, h);
                auto spec = new gltf_texture();
                spec->ldr = image4b(w, h);
                for (auto j = 0; j < h; j++) {
                    for (auto i = 0; i < w; i++) {
                        auto u = i / (float)w, v = j / (float)h;
                        auto base = vec4b{255, 255, 255, 255};
                        auto metallic = vec4b{255, 255, 255, 255};
                        if (mr->base_txt) {
                            auto ii = (int)(u * mr->base_txt->ldr.width());
                            auto jj = (int)(v * mr->base_txt->ldr.height());
                            base = mr->base_txt->ldr[{ii, jj}];
                        } else {
                            base = linear_to_srgb(vec4f(mr->base.x, mr->base.y,
                                                        mr->base.z, mr->opacity));
                        }
                        if (mr->metallic_txt) {
                            auto ii = (int)(u * mr->metallic_txt->ldr.width());
                            auto jj = (int)(v * mr->metallic_txt->ldr.height());
                            metallic = mr->metallic_txt->ldr[{ii, jj}];
                        } else {
                            metallic = linear_to_srgb(
                                vec4f(1, mr->roughness, mr->metallic, 1));
                        }
                        auto kb_txt = srgb_to_linear(base);
                        auto km_txt = srgb_to_linear(metallic);
                        auto kb = vec3f{kb_txt.x, kb_txt.y, kb_txt.z};
                        auto km = km_txt.z;
                        auto kd = kb * (1-km);
                        auto ks = kb * km + vec3f{(1 - km) * 0.04f};
                        diff->ldr[{i, j}] =
                        linear_to_srgb({kd.x, kd.y, kd.z, kb_txt.w});
                        spec->ldr[{i, j}] = linear_to_srgb(
                            {ks.x, ks.y, ks.z, 1 - km_txt.y});
                    }
                }
                txts[mr_txt] = {diff, spec};
                scns->textures.push_back(diff);
                scns->textures.push_back(spec);
            }
            std::tie(sg->diffuse_txt, sg->specular_txt) = txts[mr_txt];
        } else {
            sg->diffuse = mr->base * (1 - mr->metallic);
            sg->specular =
                mr->base * mr->metallic + vec3f{(1 - mr->metallic) * 0.04f};
            sg->opacity = mr->opacity;
            sg->glossiness = 1 - mr->roughness;
        }
    }
}

}  // namespace ygl
