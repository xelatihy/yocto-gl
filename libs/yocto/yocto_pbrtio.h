//
// # Yocto/PbrtIO: Serialization for Pbrt models
//
// Yocto/PbrtIO is a collection of utilities for loading and saving scenes
// in the Pbrt format.
// Yocto/PbrtIO is implemented in `yocto_pbrtio.h` and `yocto_pbrtio.cpp`,
// and depends on `fast_float.h` for number parsing.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef _YOCTO_PBRTIO_H_
#define _YOCTO_PBRTIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  array<float, 12> frame      = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12> frend      = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<int, 2>    resolution = {0, 0};
  float            lens       = 0;
  float            aspect     = 0;
  float            focus      = 0;
  float            aperture   = 0;
};

// Pbrt material
struct pbrt_texture {
  string          name     = "";
  array<float, 3> constant = {1, 1, 1};
  string          filename = "";
};

// Pbrt material type (simplified and only for the materials that matter here)
enum struct pbrt_mtype {
  // clang-format off
  matte, plastic, metal, glass, thinglass, subsurface
  // clang-format on
};

// Pbrt material
struct pbrt_material {
  string          name            = "";
  pbrt_mtype      type            = pbrt_mtype::matte;
  array<float, 3> emission        = {0, 0, 0};
  array<float, 3> color           = {0, 0, 0};
  float           roughness       = 0;
  float           ior             = 1.5f;
  float           opacity         = 1;
  int             color_tex       = -1;
  array<float, 3> volmeanfreepath = {0, 0, 0};
  array<float, 3> volscatter      = {0, 0, 0};
  float           volscale        = 0.01f;
};

// Pbrt shape
struct pbrt_shape {
  array<float, 12>         frame     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>         frend     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  bool                     instanced = false;
  vector<array<float, 12>> instances = {};
  vector<array<float, 12>> instaends = {};
  int                      material  = -1;
  string                   filename_ = "";
  vector<array<float, 3>>  positions = {};
  vector<array<float, 3>>  normals   = {};
  vector<array<float, 2>>  texcoords = {};
  vector<array<int, 3>>    triangles = {};
};

// Pbrt lights
struct pbrt_light {
  array<float, 12>        frame          = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>        frend          = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 3>         emission       = {0, 0, 0};
  array<float, 3>         from           = {0, 0, 0};
  array<float, 3>         to             = {0, 0, 0};
  bool                    distant        = false;
  array<float, 3>         area_emission  = {0, 0, 0};
  array<float, 12>        area_frame     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>        area_frend     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  vector<array<int, 3>>   area_triangles = {};
  vector<array<float, 3>> area_positions = {};
  vector<array<float, 3>> area_normals   = {};
};
struct pbrt_environment {
  array<float, 12> frame        = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12> frend        = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 3>  emission     = {0, 0, 0};
  int              emission_tex = -1;
};

// Pbrt model
struct pbrt_model {
  // pbrt data
  vector<string>           comments     = {};
  vector<pbrt_camera>      cameras      = {};
  vector<pbrt_shape>       shapes       = {};
  vector<pbrt_environment> environments = {};
  vector<pbrt_light>       lights       = {};
  vector<pbrt_material>    materials    = {};
  vector<pbrt_texture>     textures     = {};
};

// Load/save pbrt
[[nodiscard]] bool load_pbrt(const string& filename, pbrt_model& pbrt,
    string& error, bool ply_meshes = false);
[[nodiscard]] bool save_pbrt(const string& filename, const pbrt_model& pbrt,
    string& error, bool ply_meshes = false);

}  // namespace yocto

#endif
