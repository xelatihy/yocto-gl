//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
using namespace yocto;

#include <map>
#include <memory>
using namespace std;

#include "ext/CLI11.hpp"
#include "ext/Timer.hpp"
#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// construct a scene from io
shared_ptr<trace_scene> make_scene(shared_ptr<sceneio_model> ioscene) {
  auto scene = make_trace_scene();

  for (auto iocamera : ioscene->cameras) {
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
  }

  auto texture_map =
      unordered_map<shared_ptr<sceneio_texture>, shared_ptr<trace_texture>>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    auto texture = add_texture(scene);
    if (!iotexture->hdr.empty()) {
      set_texture(texture, std::move(iotexture->hdr));
    } else if (!iotexture->ldr.empty()) {
      set_texture(texture, std::move(iotexture->ldr));
    }
    texture_map[iotexture] = texture;
  }

  auto material_map =
      unordered_map<shared_ptr<sceneio_material>, shared_ptr<trace_material>>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    auto material = add_material(scene);
    set_emission(material, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(
        material, iomaterial->color, texture_map.at(iomaterial->color_tex));
    set_specular(material, iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_ior(material, iomaterial->ior);
    set_metallic(material, iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_transmission(material, iomaterial->transmission, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->transmission_tex));
    set_roughness(material, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(
        material, iomaterial->opacity, texture_map.at(iomaterial->opacity_tex));
    set_thin(material, iomaterial->thin);
    set_normalmap(material, texture_map.at(iomaterial->normal_tex));
    set_scattering(material, iomaterial->scattering, iomaterial->scanisotropy,
        texture_map.at(iomaterial->scattering_tex));
    material_map[iomaterial] = material;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
  }

  auto shape_map =
      unordered_map<shared_ptr<sceneio_shape>, shared_ptr<trace_shape>>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    set_quads(shape, ioshape->quads);
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_colors(shape, ioshape->colors);
    set_radius(shape, ioshape->radius);
    set_tangents(shape, ioshape->tangents);
    shape_map[ioshape] = shape;
  }

  auto instance_map =
      unordered_map<shared_ptr<sceneio_instance>, shared_ptr<trace_instance>>{};
  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    auto instance = add_instance(scene);
    set_frames(instance, ioinstance->frames);
    instance_map[ioinstance] = instance;
  }

  for (auto ioobject : ioscene->objects) {
    auto object = add_object(scene);
    set_frame(object, ioobject->frame);
    set_shape(object, shape_map.at(ioobject->shape));
    set_material(object, material_map.at(ioobject->material));
    set_instance(object, instance_map.at(ioobject->instance));
  }

  for (auto ioenvironment : ioscene->environments) {
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  return scene;
}

int run_app(int argc, const char* argv[]) {
  // options
  auto params     = trace_params{};
  auto batch      = 16;
  auto save_batch = false;
  auto add_skyenv = false;
  auto validate   = false;
  auto imfilename = "out.hdr"s;
  auto filename   = "scene.json"s;

  // maps for getting param
  auto trace_sampler_map = map<string, trace_sampler_type>{};
  for (auto idx = 0; idx < trace_sampler_names.size(); idx++) {
    trace_sampler_map[trace_sampler_names[idx]] = (trace_sampler_type)idx;
  }
  auto trace_falsecolor_map = map<string, trace_falsecolor_type>{};
  for (auto idx = 0; idx < trace_falsecolor_names.size(); idx++) {
    trace_falsecolor_map[trace_falsecolor_names[idx]] =
        (trace_falsecolor_type)idx;
  }
  auto trace_bvh_map = map<string, trace_bvh_type>{};
  for (auto idx = 0; idx < trace_bvh_names.size(); idx++) {
    trace_bvh_map[trace_bvh_names[idx]] = (trace_bvh_type)idx;
  }

  // parse command line
  auto cli = CLI::App{"Offline path tracing"};
  cli.add_option("--camera", params.camera, "Camera index.");
  cli.add_option("--resolution,-r", params.resolution, "Image resolution.");
  cli.add_option("--samples,-s", params.samples, "Number of samples.");
  cli.add_option("--tracer,-t", params.sampler, "Trace type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_map));
  cli.add_option(
         "--falsecolor,-F", params.falsecolor, "Tracer false color type.")
      ->transform(CLI::CheckedTransformer(trace_falsecolor_map));
  cli.add_option("--bounces", params.bounces, "Maximum number of bounces.");
  cli.add_option("--clamp", params.clamp, "Final pixel clamping.");
  cli.add_flag("--filter", params.tentfilter, "Filter image.");
  cli.add_option("--batch,-b", batch, "Samples per batch.");
  cli.add_flag("--env-hidden,!--no-env-hidden", params.envhidden,
      "Environments are hidden in renderer");
  cli.add_option("--save-batch", save_batch, "Save images progressively");
  cli.add_option("--bvh", params.bvh, "Bvh type")
      ->transform(CLI::CheckedTransformer(trace_bvh_map));
  cli.add_flag("--add-skyenv", add_skyenv, "Add sky envmap");
  cli.add_option("--output-image,-o", imfilename, "Image filename");
  cli.add_flag("--validate", validate, "Validate scene");
  cli.add_option("scene", filename, "Scene filename")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // scene loading
  auto ioscene = make_shared<sceneio_model>();
  {
    auto timer = CLI::AutoTimer("load");
    ioscene    = load_scene(filename);
  }

  // add components
  if (validate) {
    auto timer = CLI::AutoTimer("validate");
    for (auto& error : scene_validation(ioscene)) std::cout << error << "\n";
  }

  // convert scene
  auto scene = shared_ptr<trace_scene>();
  {
    auto timer = CLI::AutoTimer("convert");
    scene      = make_scene(ioscene);
  }

  // cleanup
  ioscene = nullptr;

  // build bvh
  {
    auto timer = CLI::AutoTimer("bvh");
    init_bvh(scene, params);
  }

  // init renderer
  {
    auto timer = CLI::AutoTimer("lights");
    init_lights(scene);
  }

  // fix renderer type if no lights
  if (scene->lights.empty() && is_sampler_lit(params)) {
    std::cout << "no lights presents, switching to eyelight shader\n";
    params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  auto state  = make_state(scene, params);
  auto render = image{state->size(), zero4f};

  // render
  for (auto sample = 0; sample < params.samples; sample += batch) {
    auto nsamples = min(batch, params.samples - sample);
    auto timer = CLI::AutoTimer{"rendering samples " + std::to_string(sample) +
                                "/" + std::to_string(params.samples)};
    render     = trace_samples(state, scene, nsamples, params);
    if (save_batch) {
      auto outfilename = fs::path(imfilename)
                             .replace_extension(
                                 "-s" + std::to_string(sample + nsamples) +
                                 fs::path(imfilename).extension().string())
                             .string();
      auto timer = CLI::AutoTimer{"saving " + outfilename};
      save_image(outfilename, render);
    }
  }

  // save image
  {
    auto timer = CLI::AutoTimer("save");
    save_image(imfilename, render);
  }

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
