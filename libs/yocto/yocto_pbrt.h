//
// # Yocto/Pbrt: Tiny library for Pbrt parsing and writing
//
// Yocto/Pbrt is a tiny library for loading and saving Pbrt.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#ifndef _YOCTO_PBRT_H_
#define _YOCTO_PBRT_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <string>
#include <vector>
#include <algorithm>
#include <memory>

// -----------------------------------------------------------------------------
// PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::pbrt {

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  frame3f frame      = identity3x4f;
  frame3f frend      = identity3x4f;
  vec2i   resolution = {0, 0};
  float   lens       = 0;
  float   aspect     = 0;
  float   focus      = 0;
  float   aperture   = 0;
};

// Pbrt material
struct pbrt_material {
  // material parameters
  std::string name            = "";
  vec3f       emission        = {0, 0, 0};
  vec3f       color           = {0, 0, 0};
  float       specular        = 0;
  float       metallic        = 0;
  float       transmission    = 0;
  float       roughness       = 0;
  float       ior             = 1.5;
  float       opacity         = 1;
  std::string color_tex       = "";
  std::string opacity_tex     = "";
  std::string alpha_tex       = "";
  bool        thin            = true;
  vec3f       volmeanfreepath = {0, 0, 0};
  vec3f       volscatter      = {0, 0, 0};
  float       volscale        = 0.01;
};

// Pbrt shape
struct pbrt_shape {
  // frames
  frame3f              frame     = identity3x4f;
  frame3f              frend     = identity3x4f;
  std::vector<frame3f> instances = {};
  std::vector<frame3f> instaends = {};
  // shape
  std::string        filename_ = "";
  std::vector<vec3f> positions = {};
  std::vector<vec3f> normals   = {};
  std::vector<vec2f> texcoords = {};
  std::vector<vec3i> triangles = {};
  // material
  pbrt::pbrt_material* material = nullptr;
};

// Pbrt lights
struct pbrt_light {
  // light parameters
  frame3f frame    = identity3x4f;
  frame3f frend    = identity3x4f;
  vec3f   emission = {0, 0, 0};
  vec3f   from     = {0, 0, 0};
  vec3f   to       = {0, 0, 0};
  bool    distant  = false;
  // arealight approximation
  vec3f              area_emission  = {0, 0, 0};
  frame3f            area_frame     = identity3x4f;
  frame3f            area_frend     = identity3x4f;
  std::vector<vec3i> area_triangles = {};
  std::vector<vec3f> area_positions = {};
  std::vector<vec3f> area_normals   = {};
};
struct pbrt_environment {
  // environment approximation
  frame3f     frame        = identity3x4f;
  frame3f     frend        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  std::string emission_tex = "";
};

// Pbrt model
struct pbrt_model {
  // pbrt data
  std::vector<std::string>        comments     = {};
  std::vector<pbrt::pbrt_camera*>      cameras      = {};
  std::vector<pbrt::pbrt_shape*>       shapes       = {};
  std::vector<pbrt::pbrt_environment*> environments = {};
  std::vector<pbrt::pbrt_light*>       lights       = {};
  std::vector<pbrt::pbrt_material*>    materials    = {};

  // cleanup
  ~pbrt_model();
};

// Load/save pbrt
bool load_pbrt(
    const std::string& filename, pbrt::pbrt_model* pbrt, std::string& error);
bool save_pbrt(const std::string& filename, pbrt::pbrt_model* pbrt,
    std::string& error, bool ply_meshes = false);

// Create pbrt
pbrt::pbrt_camera*      add_camera(pbrt::pbrt_model* pbrt);
pbrt::pbrt_shape*       add_shape(pbrt::pbrt_model* pbrt);
pbrt::pbrt_material*    add_material(pbrt::pbrt_model* pbrt);
pbrt::pbrt_environment* add_environment(pbrt::pbrt_model* pbrt);
pbrt::pbrt_light*       add_light(pbrt::pbrt_model* pbrt);

}  // namespace yocto::pbrt

#endif
