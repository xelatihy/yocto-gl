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
scene_shape* add_shape(
    scene_model* scene, const string& name, const quads_shape& shape_data) {
  auto shape       = add_shape(scene, name);
  shape->quads     = shape_data.quads;
  shape->positions = shape_data.positions;
  shape->normals   = shape_data.normals;
  shape->texcoords = shape_data.texcoords;
  return shape;
}
scene_shape* add_shape(
    scene_model* scene, const string& name, const triangles_shape& shape_data) {
  auto shape       = add_shape(scene, name);
  shape->triangles = shape_data.triangles;
  shape->positions = shape_data.positions;
  shape->normals   = shape_data.normals;
  shape->texcoords = shape_data.texcoords;
  return shape;
}
scene_shape* add_shape(
    scene_model* scene, const string& name, const lines_shape& shape_data) {
  auto shape       = add_shape(scene, name);
  shape->lines     = shape_data.lines;
  shape->positions = shape_data.positions;
  shape->normals   = shape_data.normals;
  shape->texcoords = shape_data.texcoords;
  shape->radius    = shape_data.radius;
  return shape;
}
scene_shape* add_shape(
    scene_model* scene, const string& name, const points_shape& shape_data) {
  auto shape       = add_shape(scene, name);
  shape->points    = shape_data.points;
  shape->positions = shape_data.positions;
  shape->normals   = shape_data.normals;
  shape->texcoords = shape_data.texcoords;
  shape->radius    = shape_data.radius;
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
    const vec3f& color, scene_texture* color_tex) {
  auto material       = add_material(scene, name);
  material->color     = color;
  material->color_tex = color_tex;
  return material;
}
scene_material* add_specular_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, float ior = 1.5, float specular = 1,
    scene_texture* specular_tex = nullptr, const vec3f& spectint = {1, 1, 1},
    scene_texture* spectint_tex = nullptr) {
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
  return material;
}
scene_material* add_metallic_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, float metallic = 1,
    scene_texture* metallic_tex = nullptr) {
  auto material           = add_material(scene, name);
  material->color         = color;
  material->color_tex     = color_tex;
  material->metallic      = metallic;
  material->metallic_tex  = metallic_tex;
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
  return material;
}
scene_material* add_volumetric_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, const vec3f& scattering = {0, 0, 0},
    scene_texture* scattering_tex = nullptr, float ior = 1.5,
    float specular = 1, scene_texture* specular_tex = nullptr,
    float transmission = 1, scene_texture* transmission_tex = nullptr) {
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
  material->thin             = false;
  return material;
}
scene_material* add_coated_material(scene_model* scene, const string& name,
    const vec3f& color, scene_texture* color_tex, float roughness,
    scene_texture* roughness_tex = nullptr, float ior = 1.5, float specular = 1,
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
  return material;
}
scene_material* add_bumped_material(scene_model* scene, const string& name,
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

struct test_params {
  enum struct shapes_type { features1 };
  enum struct materials_type { features1 };
  enum struct environments_type { none, sky, sunsky };
  enum struct arealights_type { none, area };
  enum struct floor_type { none, floor };
  environments_type environments = environments_type::sky;
  arealights_type   arealights   = arealights_type::area;
  floor_type        floor        = floor_type::floor;
  shapes_type       shapes       = shapes_type::features1;
  materials_type    materials    = materials_type::features1;
};

// Scene test
void make_test(scene_model* scene, const test_params& params) {
  // cameras
  add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
      {0, 1, 0}, 0.05, 2.4, 0);
  // TODO(fabio): port other cameras
  switch (params.environments) {
    case test_params::environments_type::none: break;
    case test_params::environments_type::sky:
      add_environment(scene, "sky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, false, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
      break;
    case test_params::environments_type::sunsky:
      // TODO(fabio): sunsky
      add_environment(scene, "sunsky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, true, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
      break;
  }
  switch (params.arealights) {
    case test_params::arealights_type::none: break;
    case test_params::arealights_type::area:
      add_instance(scene, "arealight1",
          lookat_frame({-0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight1", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(scene, "arealight1", {20, 20, 20}, nullptr));
      add_instance(scene, "arealight2",
          lookat_frame({+0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(scene, "arealight2", {20, 20, 20}, nullptr));
      break;
  }
  switch (params.floor) {
    case test_params::floor_type::none: break;
    case test_params::floor_type::floor:
      add_instance(scene, "floor", identity3x4f,
          add_shape(scene, "floor", make_floor({1, 1}, {2, 2}, {20, 20})),
          add_matte_material(scene, "floor", {1, 1, 1},
              add_texture(scene, "floor", make_grid({1024, 1024}))));
      break;
  }
  auto shapes    = vector<scene_shape*>{};
  auto materials = vector<scene_material*>{};
  switch (params.shapes) {
    case test_params::shapes_type::features1:
      auto bunny  = add_shape(scene, "bunny", make_sphere(32, 0.075, 1));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
      break;
  }
  switch (params.materials) {
    case test_params::materials_type::features1:
      materials = {
          // TODO(fabio): coated material
          add_coated_material(scene, "coated", {1, 1, 1},
              add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})), 0.2),
          // TODO(fabio): radius 0.2
          add_volumetric_material(scene, "glass", {1, 0.5, 0.5}, nullptr, 0),
          add_volumetric_material(scene, "jade", {0.5, 0.5, 0.5}, nullptr, 0,
              nullptr, {0.3, 0.6, 0.3}),
          add_metallic_material(
              scene, "gold", {0.66, 0.45, 0.34}, nullptr, 0.2),
          // TODO(fabio): normal "bumps-normal"
          add_bumped_material(scene, "bumped", {0.5, 0.7, 0.5}, nullptr, 0.2,
              nullptr,
              add_texture(scene, "bumps-normal",
                  bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true)),
      };
      break;
  }
  for (auto idx = 0; idx < 5; idx++) {
    add_instance(scene, shapes[idx]->name + "-" + materials[idx]->name,
        {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx - 2), 0.075, 0}},
        shapes[idx], materials[idx]);
  }
}

// Scene presets used ofr testing.
bool make_preset(scene_model* scene, const string& type, string& error) {
  if (type == "cornellbox") {
    make_cornellbox(scene);
    return true;
  } else if (type == "features1") {
    auto params      = test_params{};
    params.shapes    = test_params::shapes_type::features1;
    params.materials = test_params::materials_type::features1;
    make_test(scene, params);
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
