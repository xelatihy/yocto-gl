//
// # Yocto/CuTrace: Path tracing on Cuda/Optix
//
// Yocto/CuTrace is a simple path tracer written on the Yocto/Scene model.
// Yocto/CuTrace is implemented in `yocto_cutrace.h`, `yocto_cutrace.cpp`,
// and `yocto_cutrace.cu`.
//
// THIS IS AN EXPERIMENTAL LIBRARY THAT IS NOT READY FOR PRIME TIME
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
//

#ifndef _YOCTO_CUTRACE_H_
#define _YOCTO_CUTRACE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_scene.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Type of tracing algorithm
enum struct cutrace_sampler_type {
  path,        // path tracing
  pathdirect,  // path tracing with direct
  pathmis,     // path tracing with mis
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  eyelightao,  // eyelight with ambient occlusion
  furnace,     // furnace test
  falsecolor,  // false color rendering
};
// Type of false color visualization
enum struct cutrace_falsecolor_type {
  // clang-format off
  position, normal, frontfacing, gnormal, gfrontfacing, texcoord, mtype, color,
  emission, roughness, opacity, metallic, delta, instance, shape, material, 
  element, highlight
  // clang-format on
};

// Default trace seed
const auto cutrace_default_seed = 961748941ull;

// Options for trace functions
struct cutrace_params {
  int                     camera         = 0;
  int                     resolution     = 1280;
  cutrace_sampler_type    sampler        = cutrace_sampler_type::path;
  cutrace_falsecolor_type falsecolor     = cutrace_falsecolor_type::color;
  int                     samples        = 512;
  int                     bounces        = 8;
  float                   clamp          = 10;
  bool                    nocaustics     = false;
  bool                    envhidden      = false;
  bool                    tentfilter     = false;
  uint64_t                seed           = cutrace_default_seed;
  bool                    embreebvh      = false;
  bool                    highqualitybvh = false;
  bool                    noparallel     = false;
  int                     pratio         = 8;
  float                   exposure       = 0;
  bool                    filmic         = false;
  bool                    denoise        = false;
  int                     batch          = 1;
};

// Progressively computes an image.
image_data cutrace_image(const scene_data& scene, const cutrace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {}  // namespace yocto

// -----------------------------------------------------------------------------
// ENUM LABELS
// -----------------------------------------------------------------------------
namespace yocto {

// trace sampler names
inline const auto cutrace_sampler_names = vector<string>{"path", "pathdirect",
    "pathmis", "naive", "eyelight", "eyelightao", "furnace", "falsecolor"};

// false color names
inline const auto cutrace_falsecolor_names = vector<string>{"position",
    "normal", "frontfacing", "gnormal", "gfrontfacing", "texcoord", "mtype",
    "color", "emission", "roughness", "opacity", "metallic", "delta",
    "instance", "shape", "material", "element", "highlight"};

// trace sampler labels
inline const auto cutrace_sampler_labels =
    vector<pair<cutrace_sampler_type, string>>{
        {cutrace_sampler_type::path, "path"},
        {cutrace_sampler_type::pathdirect, "pathdirect"},
        {cutrace_sampler_type::pathmis, "pathmis"},
        {cutrace_sampler_type::naive, "naive"},
        {cutrace_sampler_type::eyelight, "eyelight"},
        {cutrace_sampler_type::eyelightao, "eyelightao"},
        {cutrace_sampler_type::furnace, "furnace"},
        {cutrace_sampler_type::falsecolor, "falsecolor"}};

// false color labels
inline const auto cutrace_falsecolor_labels =
    vector<pair<cutrace_falsecolor_type, string>>{
        {cutrace_falsecolor_type::position, "position"},
        {cutrace_falsecolor_type::normal, "normal"},
        {cutrace_falsecolor_type::frontfacing, "frontfacing"},
        {cutrace_falsecolor_type::gnormal, "gnormal"},
        {cutrace_falsecolor_type::gfrontfacing, "gfrontfacing"},
        {cutrace_falsecolor_type::texcoord, "texcoord"},
        {cutrace_falsecolor_type::mtype, "mtype"},
        {cutrace_falsecolor_type::color, "color"},
        {cutrace_falsecolor_type::emission, "emission"},
        {cutrace_falsecolor_type::roughness, "roughness"},
        {cutrace_falsecolor_type::opacity, "opacity"},
        {cutrace_falsecolor_type::metallic, "metallic"},
        {cutrace_falsecolor_type::delta, "delta"},
        {cutrace_falsecolor_type::instance, "instance"},
        {cutrace_falsecolor_type::shape, "shape"},
        {cutrace_falsecolor_type::material, "material"},
        {cutrace_falsecolor_type::element, "element"},
        {cutrace_falsecolor_type::highlight, "highlight"}};

}  // namespace yocto

#endif
