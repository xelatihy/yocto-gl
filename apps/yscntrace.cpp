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

#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
using namespace yocto;

#include <memory>
using std::make_shared;

// construct a scene from io
void init_scene(trace_scene* scene, sceneio_model* ioscene) {
  for (auto iocamera : ioscene->cameras) {
    auto camera = add_camera(scene);
    set_camera_frame(camera, iocamera->frame);
    set_camera_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_camera_focus(camera, iocamera->aperture, iocamera->focus);
  }

  auto texture_map     = unordered_map<sceneio_texture*, trace_texture*>{};
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

  auto material_map     = unordered_map<sceneio_material*, trace_material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    auto material = add_material(scene);
    set_shape_emission(material, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_shape_color(
        material, iomaterial->color, texture_map.at(iomaterial->color_tex));
    set_shape_specular(material, iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_shape_ior(material, iomaterial->ior);
    set_shape_metallic(material, iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_shape_transmission(material, iomaterial->transmission, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->transmission_tex));
    set_shape_roughness(material, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_shape_opacity(
        material, iomaterial->opacity, texture_map.at(iomaterial->opacity_tex));
    set_shape_thin(material, iomaterial->thin);
    set_shape_normalmap(material, texture_map.at(iomaterial->normal_tex));
    set_shape_scattering(material, iomaterial->scattering,
        iomaterial->scanisotropy, texture_map.at(iomaterial->scattering_tex));
    material_map[iomaterial] = material;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
  }

  auto shape_map     = unordered_map<sceneio_shape*, trace_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    auto shape = add_shape(scene);
    set_shape_points(shape, ioshape->points);
    set_shape_lines(shape, ioshape->lines);
    set_shape_triangles(shape, ioshape->triangles);
    set_shape_quads(shape, ioshape->quads);
    set_shape_positions(shape, ioshape->positions);
    set_shape_normals(shape, ioshape->normals);
    set_shape_texcoords(shape, ioshape->texcoords);
    set_shape_colors(shape, ioshape->colors);
    set_shape_radius(shape, ioshape->radius);
    set_shape_tangents(shape, ioshape->tangents);
    shape_map[ioshape] = shape;
  }

  auto instance_map     = unordered_map<sceneio_instance*, trace_instance*>{};
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
    set_environment_frame(environment, ioenvironment->frame);
    set_environment_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }
}

void run_app(int argc, const char* argv[]) {
  // options
  auto params     = trace_params{};
  auto batch      = 16;
  auto save_batch = false;
  auto add_skyenv = false;
  auto validate   = false;
  auto imfilename = "out.hdr"s;
  auto filename   = "scene.json"s;

  // parse command line
  auto cli = make_cli("yscntrace", "Offline path tracing");
  add_cli_option(cli, "--camera", params.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", params.resolution, "Image resolution.");
  add_cli_option(cli, "--samples,-s", params.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)params.sampler, "Trace type.",
      trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", params.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", params.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", params.tentfilter, "Filter image.");
  add_cli_option(cli, "--batch,-b", batch, "Samples per batch.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", params.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(cli, "--save-batch", save_batch, "Save images progressively");
  add_cli_option(cli, "--bvh", (int&)params.bvh, "Bvh type", trace_bvh_names);
  add_cli_option(cli, "--add-skyenv", add_skyenv, "Add sky envmap");
  add_cli_option(cli, "--output-image,-o", imfilename, "Image filename");
  add_cli_option(cli, "--validate", validate, "Validate scene");
  add_cli_option(cli, "scene", filename, "Scene filename", true);
  parse_cli(cli, argc, argv);

  // scene loading
  auto ioscene    = make_shared<sceneio_model>();
  auto load_timer = print_timed("loading scene");
  load_scene(filename, ioscene.get());
  print_elapsed(load_timer);

  // add components
  if (validate) {
    auto validate_timer = print_timed("validating");
    auto errors         = scene_validation(ioscene.get());
    for (auto& error : errors) print_info(error);
    print_elapsed(validate_timer);
  }

  // convert scene
  auto convert_timer = print_timed("converting");
  auto scene         = make_shared<trace_scene>();
  init_scene(scene.get(), ioscene.get());
  print_elapsed(convert_timer);

  // cleanup
  ioscene = nullptr;

  // build bvh
  auto bvh_timer = print_timed("building bvh");
  init_bvh(scene.get(), params);
  print_elapsed(bvh_timer);

  // init renderer
  auto lights_timer = print_timed("building lights");
  init_lights(scene.get());
  print_elapsed(lights_timer);

  // fix renderer type if no lights
  if (scene->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, switching to eyelight shader");
    params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  auto state = make_shared<trace_state>();
  init_state(state.get(), scene.get(), params);
  auto render = image{state->size(), zero4f};

  // render
  for (auto sample = 0; sample < params.samples; sample += batch) {
    auto nsamples    = min(batch, params.samples - sample);
    auto batch_timer = print_timed("rendering samples " +
                                   std::to_string(sample) + "/" +
                                   std::to_string(params.samples));
    render = trace_samples(state.get(), scene.get(), nsamples, params);
    print_elapsed(batch_timer);
    if (save_batch) {
      auto outfilename = replace_extension(imfilename,
          "-s" + std::to_string(sample + nsamples) + get_extension(imfilename));
      save_image(outfilename, render);
    }
  }

  // save image
  auto save_timer = print_timed("saving image");
  save_image(imfilename, render);
  print_elapsed(save_timer);
}

int main(int argc, const char* argv[]) {
  try {
    run_app(argc, argv);
    return 0;
  } catch (std::exception& e) {
    print_fatal(e.what());
    return 1;
  }
}
