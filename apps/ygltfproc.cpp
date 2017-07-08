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

#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_utils.h"

#include "../yocto/ext/json.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>

using namespace nlohmann;

int main(int argc, char** argv) {
    // command line params
    auto parser =
        yu::cmdline::make_parser(argc, argv, "prints info about gltf");
    auto use_scene = parse_flag(
        parser, "--use-flag", "-s", "convert to/from scene internally");
    auto save_info =
        parse_opts(parser, "--info-filename", "-i", "info filename (JSON)", "");
    auto save_out =
        parse_opts(parser, "--out-filename", "-o", "output filename", "");
    auto filename = parse_args(parser, "filename", "input filename", "", true);
    check_parser(parser);

    // load gltf
    auto gltf = ygltf::load_gltf(filename, true);

    // scene
    if (use_scene) {
        auto scns = ygltf::gltf_to_scenes(gltf);
        delete gltf;
        auto buffer_uri = yu::path::get_basename(filename) + ".bin";
        gltf = ygltf::scenes_to_gltf(scns, buffer_uri);
        delete scns;
    }

    // saving output
    if (!save_out.empty()) ygltf::save_gltf(save_out, gltf, true, false);

    // info
    if (!save_info.empty()) {
        // start statistics
        auto info = json::object();

        // number of elements
        auto objects = json::object();
        objects["accessors"] = gltf->accessors.size();
        objects["animations"] = gltf->animations.size();
        objects["buffers"] = gltf->buffers.size();
        objects["bufferViews"] = gltf->bufferViews.size();
        objects["cameras"] = gltf->cameras.size();
        objects["images"] = gltf->images.size();
        objects["materials"] = gltf->materials.size();
        objects["meshes"] = gltf->meshes.size();
        objects["nodes"] = gltf->nodes.size();
        objects["samplers"] = gltf->samplers.size();
        objects["scenes"] = gltf->scenes.size();
        objects["skins"] = gltf->skins.size();
        objects["textures"] = gltf->textures.size();
        info["objects"] = objects;

        static const auto type_map =
            std::map<ygltf::glTFAccessorType, std::string>{
                {ygltf::glTFAccessorType::NotSet, "NotSet"},
                {ygltf::glTFAccessorType::Scalar, "Scalar"},
                {ygltf::glTFAccessorType::Vec2, "Vec2"},
                {ygltf::glTFAccessorType::Vec3, "Vec3"},
                {ygltf::glTFAccessorType::Vec4, "Vec4"},
                {ygltf::glTFAccessorType::Mat2, "Mat2"},
                {ygltf::glTFAccessorType::Mat3, "Mat3"},
                {ygltf::glTFAccessorType::Mat4, "Mat4"},
            };

        static const auto ctype_map =
            std::map<ygltf::glTFAccessorComponentType, std::string>{
                {ygltf::glTFAccessorComponentType::NotSet, "NotSet"},
                {ygltf::glTFAccessorComponentType::Byte, "Byte"},
                {ygltf::glTFAccessorComponentType::UnsignedByte,
                    "Unsigned Byte"},
                {ygltf::glTFAccessorComponentType::Short, "Short"},
                {ygltf::glTFAccessorComponentType::UnsignedShort,
                    "Unsigned Short"},
                {ygltf::glTFAccessorComponentType::UnsignedInt, "UnsignedInt"},
                {ygltf::glTFAccessorComponentType::Float, "Float"}};

        // buffer descriptions
        auto descrs = json::array();
        for (auto gdescr : ygltf::gen_buffer_descriptors(gltf)) {
            auto descr = json::object();
            descr["name"] = gdescr->name;
            descr["uri"] = gdescr->uri;
            descr["size"] = gdescr->size;
            descr["sections"] = json::array();
            for (auto gsect : gdescr->sections) {
                auto sect = json::object();
                sect["refcount"] = gsect->refcount;
                sect["start"] = gsect->start;
                sect["size"] = gsect->size;
                sect["stride"] = gsect->stride;
                sect["count"] = gsect->count;
                sect["type"] = type_map.at(gsect->type);
                sect["ctype"] = ctype_map.at(gsect->ctype);
                sect["ncomp"] = gsect->ncomp;
                sect["csize"] = gsect->csize;
                descr["sections"].push_back(sect);
            }
            descrs.push_back(descr);
        }
        info["descriptors"] = descrs;

        // dump json
        auto info_txt = info.dump(4);
        yu::file::save_txtfile(save_info.c_str(), info_txt);
    }

    // cleanup
    delete gltf;
    return 0;
}
