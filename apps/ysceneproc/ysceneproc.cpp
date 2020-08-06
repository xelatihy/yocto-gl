//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

#include <memory>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

image<vec3b> rgba_to_rgb(const image<vec4b>& rgba) {
  auto rgb = image<vec3b>{rgba.imsize()};
  for (auto i = 0; i < rgb.count(); i++)
    rgb[i] = {rgba[i].x, rgba[i].y, rgba[i].z};
  return rgb;
}
image<vec3f> rgba_to_rgb(const image<vec4f>& rgba) {
  auto rgb = image<vec3f>{rgba.imsize()};
  for (auto i = 0; i < rgb.count(); i++)
    rgb[i] = {rgba[i].x, rgba[i].y, rgba[i].z};
  return rgb;
}
image<byte> rgba_to_r(const image<vec4b>& rgba) {
  auto r = image<byte>{rgba.imsize()};
  for (auto i = 0; i < r.count(); i++) r[i] = rgba[i].x;
  return r;
}
image<float> rgba_to_r(const image<vec4f>& rgba) {
  auto r = image<float>{rgba.imsize()};
  for (auto i = 0; i < r.count(); i++) r[i] = rgba[i].x;
  return r;
}

#include "yshapedata.h"
triangles_shape make_bunny(float scale = 1, bool align_middle = true) {
  auto shape      = triangles_shape{};
  shape.triangles = bunny_triangles;
  shape.positions = bunny_positions;
  shape.normals   = bunny_normals;
  shape.texcoords = bunny_texcoords;
  // scale to height 1
  auto bbox = invalidb3f;
  for (auto& t : shape.triangles) {
    bbox = merge(bbox, triangle_bounds(shape.positions[t.x],
                           shape.positions[t.y], shape.positions[t.z]));
  }
  auto yscale = 2 / size(bbox).y;
  for (auto& p : shape.positions) p *= yscale;
  if (align_middle) {
    for (auto& p : shape.positions) p.y -= 1;
  }
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}

scene_camera* add_camera(scene_model* scene, const string& name,
    const vec3f& from, const vec3f& to, const vec3f& up, float lens,
    float aspect, float aperture = 0, bool ortho = false, float film = 0.036) {
  auto camera = add_camera(scene, name);
  set_frame(camera, lookat_frame(from, to, up));
  set_lens(camera, lens, aspect, film, ortho);
  set_focus(camera, aperture, length(from - to));
  return camera;
}
scene_camera* add_camera(scene_model* scene, const string& name,
    const frame3f& frame, float lens, float aspect, float aperture = 0,
    float focus = 10, bool ortho = false, float film = 0.036) {
  auto camera = add_camera(scene, name);
  set_frame(camera, frame);
  set_lens(camera, lens, aspect, film, ortho);
  set_focus(camera, aperture, focus);
  return camera;
}
scene_instance* add_instance(scene_model* scene, const string& name,
    const frame3f& frame, scene_shape* shape, scene_material* material) {
  auto instance = add_instance(scene, name);
  set_frame(instance, frame);
  set_shape(instance, shape);
  set_material(instance, material);
  return instance;
}
scene_environment* add_environment(scene_model* scene, const string& name,
    const frame3f& frame, const vec3f& emission,
    scene_texture* emission_tex = nullptr) {
  auto environment = add_environment(scene, name);
  set_frame(environment, frame);
  set_emission(environment, emission, emission_tex);
  return environment;
}
scene_texture* add_texture(scene_model* scene, const string& name,
    const image<vec4f>& img, bool hdr = false, bool ldr_linear = false,
    bool single_channel = false) {
  auto texture = add_texture(scene, name);
  if (hdr) {
    if (single_channel) {
      set_texture(texture, rgba_to_r(img));
    } else {
      set_texture(texture, rgba_to_rgb(img));
    }
  } else {
    auto imgb = ldr_linear ? float_to_byte(img) : rgb_to_srgbb(img);
    if (single_channel) {
      set_texture(texture, rgba_to_r(imgb));
    } else {
      set_texture(texture, rgba_to_rgb(imgb));
    }
  }
  return texture;
}
scene_shape* add_shape(scene_model* scene, const string& name,
    const quads_shape& shape_data, int subdivisions = 0, float displacement = 0,
    scene_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->quads            = shape_data.quads;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
scene_shape* add_shape(scene_model* scene, const string& name,
    const quads_fvshape& shape_data, int subdivisions = 0,
    float displacement = 0, scene_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->quadspos         = shape_data.quadspos;
  shape->quadsnorm        = shape_data.quadsnorm;
  shape->quadstexcoord    = shape_data.quadstexcoord;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
scene_shape* add_shape(scene_model* scene, const string& name,
    const triangles_shape& shape_data, int subdivisions = 0,
    float displacement = 0, scene_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->triangles        = shape_data.triangles;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
scene_shape* add_shape(scene_model* scene, const string& name,
    const lines_shape& shape_data, int subdivisions = 0, float displacement = 0,
    scene_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->lines            = shape_data.lines;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->radius           = shape_data.radius;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
scene_shape* add_shape(scene_model* scene, const string& name,
    const points_shape& shape_data, int subdivisions = 0,
    float displacement = 0, scene_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->points           = shape_data.points;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->radius           = shape_data.radius;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
scene_material* add_emission_material(scene_model* scene, const string& name,
    const vec3f& emission, scene_texture* emission_tex) {
  auto material          = add_material(scene, name);
  material->emission     = emission;
  material->emission_tex = emission_tex;
  return material;
}
scene_material* add_matte_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex,
    scene_texture* normal_tex = nullptr) {
  auto material        = add_material(scene, name);
  material->color      = color;
  material->color_tex  = color_tex;
  material->roughness  = 1;
  material->normal_tex = normal_tex;
  return material;
}
scene_material* add_specular_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, scene_texture* normal_tex = nullptr,
    float ior = 1.5, float specular = 1, scene_texture* specular_tex = nullptr,
    const vec3f& spectint = {1, 1, 1}, scene_texture* spectint_tex = nullptr) {
  auto material           = add_material(scene, name);
  material->color         = color;
  material->color_tex     = color_tex;
  material->specular      = specular;
  material->specular_tex  = specular_tex;
  material->spectint      = spectint;
  material->spectint_tex  = spectint_tex;
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
  material->ior           = ior;
  material->normal_tex    = normal_tex;
  return material;
}
scene_material* add_metallic_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, scene_texture* normal_tex = nullptr,
    float metallic = 1, scene_texture* metallic_tex = nullptr) {
  auto material           = add_material(scene, name);
  material->color         = color;
  material->color_tex     = color_tex;
  material->metallic      = metallic;
  material->metallic_tex  = metallic_tex;
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
  material->normal_tex    = normal_tex;
  return material;
}
scene_material* add_transmission_material(scene_model* scene,
    const string& name, const vec3f& color, scene_texture* color_tex,
    float roughness, scene_texture* roughness_tex = nullptr,
    scene_texture* normal_tex = nullptr, float ior = 1.5, float specular = 1,
    scene_texture* specular_tex = nullptr, float transmission = 1,
    scene_texture* transmission_tex = nullptr) {
  auto material              = add_material(scene, name);
  material->color            = color;
  material->color_tex        = color_tex;
  material->specular         = specular;
  material->specular_tex     = specular_tex;
  material->transmission     = transmission;
  material->transmission_tex = transmission_tex;
  material->roughness        = roughness;
  material->roughness_tex    = roughness_tex;
  material->ior              = ior;
  material->thin             = true;
  material->normal_tex       = normal_tex;
  return material;
}
scene_material* add_volumetric_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, const vec3f& scattering = {0, 0, 0},
    scene_texture* scattering_tex = nullptr,
    scene_texture* normal_tex = nullptr, float ior = 1.5,
    float scanisotropy = 0, float trdepth = 0.01, float specular = 1,
    scene_texture* specular_tex = nullptr, float transmission = 1,
    scene_texture* transmission_tex = nullptr) {
  auto material              = add_material(scene, name);
  material->color            = color;
  material->color_tex        = color_tex;
  material->specular         = specular;
  material->specular_tex     = specular_tex;
  material->transmission     = transmission;
  material->transmission_tex = transmission_tex;
  material->roughness        = roughness;
  material->roughness_tex    = roughness_tex;
  material->scattering       = scattering;
  material->scattering_tex   = scattering_tex;
  material->ior              = ior;
  material->scanisotropy     = scanisotropy;
  material->trdepth          = trdepth;
  material->normal_tex       = normal_tex;
  material->thin             = false;
  return material;
}
scene_material* add_volumetrict_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, const vec3f& scattering = {0, 0, 0},
    scene_texture* scattering_tex = nullptr,
    scene_texture* normal_tex = nullptr, float ior = 1.5,
    float scanisotropy = 0, float trdepth = 0.01, float specular = 1,
    scene_texture* specular_tex = nullptr, float translucency = 1,
    scene_texture* translucency_tex = nullptr) {
  auto material              = add_material(scene, name);
  material->color            = color;
  material->color_tex        = color_tex;
  material->specular         = specular;
  material->specular_tex     = specular_tex;
  material->translucency     = translucency;
  material->translucency_tex = translucency_tex;
  material->roughness        = roughness;
  material->roughness_tex    = roughness_tex;
  material->scattering       = scattering;
  material->scattering_tex   = scattering_tex;
  material->ior              = ior;
  material->scanisotropy     = scanisotropy;
  material->trdepth          = trdepth;
  material->normal_tex       = normal_tex;
  material->thin             = false;
  return material;
}
scene_material* add_specular_coated_material(scene_model* scene,
    const string& name, const vec3f& color, scene_texture* color_tex,
    float roughness, scene_texture* roughness_tex = nullptr,
    scene_texture* normal_tex = nullptr, float ior = 1.5, float specular = 1,
    scene_texture* specular_tex = nullptr, float coat = 1,
    scene_texture* coat_tex = nullptr) {
  auto material           = add_material(scene, name);
  material->color         = color;
  material->color_tex     = color_tex;
  material->specular      = specular;
  material->specular_tex  = specular_tex;
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
  material->coat          = coat;
  material->coat_tex      = coat_tex;
  material->ior           = ior;
  material->normal_tex    = normal_tex;
  return material;
}
scene_material* add_metallic_coated_material(scene_model* scene,
    const string& name, const vec3f& color, scene_texture* color_tex,
    float roughness, scene_texture* roughness_tex = nullptr,
    scene_texture* normal_tex = nullptr, float metallic = 1,
    scene_texture* metallic_tex = nullptr, float coat = 1,
    scene_texture* coat_tex = nullptr) {
  auto material           = add_material(scene, name);
  material->color         = color;
  material->color_tex     = color_tex;
  material->metallic      = metallic;
  material->metallic_tex  = metallic_tex;
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
  material->coat          = coat;
  material->coat_tex      = coat_tex;
  material->normal_tex    = normal_tex;
  return material;
}
scene_material* add_transparent_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float opacity = 1,
    scene_texture* normal_tex = nullptr) {
  auto material        = add_material(scene, name);
  material->color      = color;
  material->color_tex  = color_tex;
  material->roughness  = 1;
  material->opacity    = opacity;
  material->normal_tex = normal_tex;
  return material;
}

enum struct test_cameras_type { standard, wide };
enum struct test_environments_type { none, sky, sunsky };
enum struct test_arealights_type { none, standard, large };
enum struct test_floor_type { none, standard };
enum struct test_instance_name_type { material, shape };
enum struct test_shapes_type {
  // clang-format off
  features1, features2, rows, bunny_sphere,
  shapes1, shapes2, shapes3
  // clang-format off
};
enum struct test_materials_type {
  // clang-format off
  features1, features2, uvgrid, hair, plastic_metal,
  materials1, materials2, materials3, materials4, materials5,
  // clang-format on
};

struct test_params {
  test_cameras_type       cameras       = test_cameras_type::standard;
  test_environments_type  environments  = test_environments_type::sky;
  test_arealights_type    arealights    = test_arealights_type::standard;
  test_floor_type         floor         = test_floor_type::standard;
  test_shapes_type        shapes        = test_shapes_type::features1;
  test_materials_type     materials     = test_materials_type::features1;
  test_instance_name_type instance_name = test_instance_name_type::material;
};

// Scene test
void make_test(scene_model* scene, const test_params& params) {
  // cameras
  switch (params.cameras) {
    case test_cameras_type::standard: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
    // TODO(fabio): fix wide camera
    case test_cameras_type::wide: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
  }
  // TODO(fabio): port other cameras
  switch (params.environments) {
    case test_environments_type::none: break;
    case test_environments_type::sky: {
      add_environment(scene, "sky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, false, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
    } break;
    case test_environments_type::sunsky: {
      add_environment(scene, "sunsky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, true, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
    } break;
  }
  switch (params.arealights) {
    case test_arealights_type::none: break;
    case test_arealights_type::standard: {
      add_instance(scene, "arealight1",
          lookat_frame({-0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight1", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(scene, "arealight1", {20, 20, 20}, nullptr));
      add_instance(scene, "arealight2",
          lookat_frame({+0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(scene, "arealight2", {20, 20, 20}, nullptr));
    } break;
    case test_arealights_type::large: {
      add_instance(scene, "largearealight1",
          lookat_frame({-0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "largearealight1", make_rect({1, 1}, {0.4, 0.4})),
          add_emission_material(
              scene, "largearealight1", {10, 10, 10}, nullptr));
      add_instance(scene, "largearealight2",
          lookat_frame({+0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "largearealight2", make_rect({1, 1}, {0.4, 0.4})),
          add_emission_material(
              scene, "largearealight2", {10, 10, 10}, nullptr));
    } break;
  }
  switch (params.floor) {
    case test_floor_type::none: break;
    case test_floor_type::standard: {
      add_instance(scene, "floor", identity3x4f,
          add_shape(scene, "floor", make_floor({1, 1}, {2, 2}, {20, 20})),
          add_matte_material(scene, "floor", {1, 1, 1},
              add_texture(scene, "floor", make_grid({1024, 1024}))));
    } break;
  }
  auto shapes = vector<scene_shape*>{}, shapesi = vector<scene_shape*>{};
  auto materials = vector<scene_material*>{};
  switch (params.shapes) {
    case test_shapes_type::features1: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::features2: {
      shapes  = {add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
          add_shape(scene, "suzanne", make_monkey(0.075f * 0.8f), 2),
          add_shape(scene, "hair",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
          add_shape(scene, "displaced", make_sphere(128, 0.075f, 1), 0, 0.025,
              add_texture(scene, "bumps-displacement", make_bumps({1024, 1024}),
                  false, true)),
          add_shape(scene, "cube",
              make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1, 1},
                  0.3 * 0.075f))};
      shapesi = {nullptr, nullptr,
          add_shape(scene, "hairi", make_sphere(32, 0.075f * 0.8f, 1)), nullptr,
          nullptr};
    } break;
    case test_shapes_type::rows: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, bunny, bunny, bunny, bunny, sphere, sphere, sphere,
          sphere, sphere};
    } break;
    case test_shapes_type::bunny_sphere: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::shapes1: {
      shapes = {
          add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
          add_shape(scene, "uvsphere-flipcap",
              make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075)),
          add_shape(scene, "disk", make_disk(32, 0.075f, 1)),
          add_shape(scene, "uvcylinder",
              make_rounded_uvcylinder(
                  {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075)),
          add_shape(scene, "cube",
              make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1, 1},
                  0.3 * 0.075f)),
      };
    } break;
    case test_shapes_type::shapes2: {
      shapes = {
          add_shape(scene, "cube-subdiv", make_fvcube(0.075), 4),
          add_shape(scene, "suzanne-subdiv", make_monkey(0.075), 2),
          add_shape(scene, "displaced", make_sphere(128, 0.075f, 1), 0, 0.025,
              add_texture(scene, "bumps-displacement", make_bumps({1024, 1024}),
                  false, true)),
          add_shape(scene, "bunny", make_bunny(0.075)),
          add_shape(scene, "teapot", make_sphere(32, 0.075, 1)),
      };
    } break;
    case test_shapes_type::shapes3: {
      shapes = {
          nullptr,
          add_shape(scene, "hair1",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
          add_shape(scene, "hair2",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f})),
          add_shape(scene, "hair3",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128})),
          nullptr,
      };
    } break;
  }
  switch (params.materials) {
    case test_materials_type::features1: {
      materials = {
          add_specular_coated_material(scene, "coated", {1, 1, 1},
              add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})), 0.2),
          add_volumetric_material(scene, "glass", {1, 0.5, 0.5}, nullptr, 0),
          add_volumetric_material(scene, "jade", {0.5, 0.5, 0.5}, nullptr, 0,
              nullptr, {0.3, 0.6, 0.3}),
          add_specular_material(scene, "bumped", {0.5, 0.7, 0.5}, nullptr, 0.2,
              nullptr,
              add_texture(scene, "bumps-normal",
                  bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true)),
          add_metallic_material(
              scene, "metal", {0.66, 0.45, 0.34}, nullptr, 0.2),
      };
    } break;
    case test_materials_type::features2: {
      auto uvgrid  = add_specular_material(scene, "uvgrid", {1, 1, 1},
          add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})), 0.2);
      auto plastic = add_specular_material(
          scene, "plastic", {0.5, 0.7, 0.5}, nullptr, 0.2);
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7}, nullptr);
      materials = {uvgrid, plastic, hair, plastic, uvgrid};
    } break;
    case test_materials_type::uvgrid: {
      auto uvgrid = add_specular_material(scene, "uvgrid", {1, 1, 1},
          add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})), 0.2);
      materials   = {uvgrid, uvgrid, uvgrid, uvgrid, uvgrid};
    } break;
    case test_materials_type::hair: {
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7}, nullptr);
      materials = {hair, hair, hair, hair, hair};
    } break;
    case test_materials_type::plastic_metal: {
      materials = {
          add_specular_material(
              scene, "plastic1", {0.5, 0.5, 0.7}, nullptr, 0.01),
          add_specular_material(
              scene, "plastic2", {0.5, 0.7, 0.5}, nullptr, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}, nullptr),
          add_metallic_material(scene, "metal1", {0.7, 0.7, 0.7}, nullptr, 0),
          add_metallic_material(
              scene, "metal2", {0.66, 0.45, 0.34}, nullptr, 0.2),
      };
    } break;
    case test_materials_type::materials1: {
      materials = {
          add_specular_material(
              scene, "plastic1", {0.5, 0.5, 0.7}, nullptr, 0.01),
          add_specular_material(
              scene, "plastic2", {0.5, 0.7, 0.5}, nullptr, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}, nullptr),
          add_metallic_material(scene, "metal1", {0.7, 0.7, 0.7}, nullptr, 0),
          add_metallic_material(
              scene, "metal2", {0.66, 0.45, 0.34}, nullptr, 0.2),
      };
    } break;
    case test_materials_type::materials2: {
      materials = {
          add_volumetric_material(scene, "glass1", {1, 1, 1}, nullptr, 0),
          add_volumetric_material(scene, "glass2", {1, 0.7, 0.7}, nullptr, 0.1),
          add_transparent_material(
              scene, "transparent", {0.7, 0.5, 0.5}, nullptr, 0.2),
          add_transmission_material(scene, "tglass1", {1, 1, 1}, nullptr, 0),
          add_transmission_material(
              scene, "tglass2", {1, 0.7, 0.7}, nullptr, 0.1),
      };
    } break;
    case test_materials_type::materials3: {
      auto bumps_normal = add_texture(scene, "bumps-normal",
          bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true);
      materials         = {
          add_specular_material(scene, "plastic1", {0.5, 0.5, 0.7}, nullptr,
              0.01, nullptr, bumps_normal),
          add_specular_coated_material(
              scene, "plastic2", {0.5, 0.7, 0.5}, nullptr, 0.2),
          add_metallic_material(scene, "metal1", {0.7, 0.7, 0.7}, nullptr, 0,
              nullptr, bumps_normal),
          add_metallic_coated_material(
              scene, "metal2", {0.66, 0.45, 0.34}, nullptr, 0.2),
          add_metallic_material(
              scene, "metal3", {0.66, 0.45, 0.34}, nullptr, 0.2),
      };
    } break;
    case test_materials_type::materials4: {
      materials = {
          add_volumetric_material(scene, "cloud", {0.65, 0.65, 0.65}, nullptr,
              0, nullptr, {0.9, 0.9, 0.9}, nullptr, nullptr, 1),
          add_volumetric_material(scene, "glass", {1, 0.5, 0.5}, nullptr, 0),
          add_volumetric_material(scene, "jade", {0.5, 0.5, 0.5}, nullptr, 0,
              nullptr, {0.3, 0.6, 0.3}),
          add_volumetrict_material(scene, "jade2", {0.5, 0.5, 0.5}, nullptr, 0,
              nullptr, {0.3, 0.6, 0.3}),
          add_volumetric_material(scene, "smoke", {0.5, 0.5, 0.5}, nullptr, 0.2,
              nullptr, {0.2, 0.2, 0.2}),
      };
    } break;
    case test_materials_type::materials5: {
      materials = {
          add_volumetric_material(scene, "skin1a", {0.76, 0.48, 0.23}, nullptr,
              0.25, nullptr, {0.436, 0.227, 0.131}, nullptr, nullptr, 1.5, -0.8,
              0.001),
          add_volumetric_material(scene, "skin2a", {0.82, 0.55, 0.4}, nullptr,
              0.25, nullptr, {0.623, 0.433, 0.343}, nullptr, nullptr, 1.5, -0.8,
              0.001),
          add_volumetric_material(scene, "skins", {0.76, 0.48, 0.23}, nullptr,
              0, nullptr, {0.436, 0.227, 0.131}, nullptr, nullptr, 1.5, -0.8,
              0.001),
          add_volumetrict_material(scene, "skin1b", {0.76, 0.48, 0.23}, nullptr,
              0.25, nullptr, {0.436, 0.227, 0.131}, nullptr, nullptr, 1.5, -0.8,
              0.001),
          add_volumetrict_material(scene, "skin2b", {0.82, 0.55, 0.4}, nullptr,
              0.25, nullptr, {0.623, 0.433, 0.343}, nullptr, nullptr, 1.5, -0.8,
              0.001),
      };
    } break;
  }
  for (auto idx = 0; idx < shapes.size(); idx++) {
    if (!shapes[idx]) continue;
    if (shapes.size() > 5) {
      add_instance(scene, shapes[idx]->name + "-" + materials[idx % 5]->name,
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
              {0.2f * (idx % 5 - 2), 0.075, -0.4f * (idx / 5)}},
          shapes[idx], materials[idx % 5]);
    } else {
      auto name = params.instance_name == test_instance_name_type::material
                      ? materials[idx]->name
                      : shapes[idx]->name;
      add_instance(scene, name,
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx % 5 - 2), 0.075, 0}},
          shapes[idx], materials[idx]);
    }
    if (!shapesi.empty() && shapesi[idx]) {
      add_instance(scene, shapesi[idx]->name,
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx - 2), 0.075, 0}},
          shapesi[idx], materials[idx]);
    }
  }
}

// Scene presets used ofr testing.
bool make_preset(scene_model* scene, const string& type, string& error) {
  if (type == "cornellbox") {
    make_cornellbox(scene);
    return true;
  } else if (type == "features1") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features1, test_materials_type::features1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "features2") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features2, test_materials_type::features2,
                   test_instance_name_type::shape});
    return true;
  } else if (type == "materials1") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials2") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials2,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials3") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials3,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials4") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials4,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials5") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials5,
                   test_instance_name_type::material});
    return true;
  } else if (type == "shapes1") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes1, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes2") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes2, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes3") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes3, test_materials_type::hair,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "environments1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "environments2") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sunsky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "arealights1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::none,
            test_arealights_type::standard, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

void make_dir(const string& dirname) {
  if (sfs::exists(dirname)) return;
  try {
    sfs::create_directories(dirname);
  } catch (...) {
    print_fatal("cannot create directory " + dirname);
  }
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto validate  = false;
  auto info      = false;
  auto copyright = ""s;
  auto output    = "out.json"s;
  auto filename  = "scene.json"s;

  // parse command line
  auto cli = make_cli("yscnproc", "Process scene");
  add_option(cli, "--info,-i", info, "print scene info");
  add_option(cli, "--copyright,-c", copyright, "copyright string");
  add_option(cli, "--validate/--no-validate", validate, "Validate scene");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "scene", filename, "input scene", true);
  parse_cli(cli, argc, argv);

  // load scene
  auto ext         = sfs::path(filename).extension().string();
  auto basename    = sfs::path(filename).stem().string();
  auto scene_guard = std::make_unique<scene_model>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (ext == ".ypreset") {
    print_progress("make preset", 0, 1);
    if (!make_preset(scene, basename, ioerror)) print_fatal(ioerror);
    print_progress("make preset", 1, 1);
  } else {
    if (!load_scene(filename, scene, ioerror, print_progress))
      print_fatal(ioerror);
  }

  // copyright
  if (copyright != "") {
    scene->copyright = copyright;
  }

  // validate scene
  if (validate) {
    for (auto& error : scene_validation(scene)) print_info("error: " + error);
  }

  // print info
  if (info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // tesselate if needed
  if (sfs::path(output).extension() != ".json") {
    tesselate_shapes(scene, print_progress);
  }

  // make a directory if needed
  make_dir(sfs::path(output).parent_path());
  if (!scene->shapes.empty())
    make_dir(sfs::path(output).parent_path() / "shapes");
  if (!scene->textures.empty())
    make_dir(sfs::path(output).parent_path() / "textures");

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // done
  return 0;
}
