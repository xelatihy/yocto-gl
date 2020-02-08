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
using std::make_unique;
using std::unique_ptr;

// construct a scene from io
void init_scene(trace_scene& scene, sceneio_model& ioscene) {
  scene = trace_scene{};

  for (auto& iocamera : ioscene.cameras) {
    auto id = add_camera(scene);
    set_camera_frame(scene, id, iocamera.frame);
    set_camera_lens(scene, id, iocamera.lens, iocamera.aspect, iocamera.film);
    set_camera_focus(scene, id, iocamera.aperture, iocamera.focus);
  }

  for (auto& iotexture : ioscene.textures) {
    auto id = add_texture(scene);
    if (!iotexture.hdr.empty()) {
      set_texture(scene, id, std::move(iotexture.hdr));
    } else if (!iotexture.ldr.empty()) {
      set_texture(scene, id, std::move(iotexture.ldr));
    }
  }

  for (auto& iosubdiv : ioscene.subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
    iosubdiv = {};
  }

  for (auto& ioshape : ioscene.shapes) {
    auto id = add_shape(scene);
    set_shape_points(scene, id, ioshape.points);
    set_shape_lines(scene, id, ioshape.lines);
    set_shape_triangles(scene, id, ioshape.triangles);
    set_shape_quads(scene, id, ioshape.quads);
    set_shape_positions(scene, id, ioshape.positions);
    set_shape_normals(scene, id, ioshape.normals);
    set_shape_texcoords(scene, id, ioshape.texcoords);
    set_shape_colors(scene, id, ioshape.colors);
    set_shape_radius(scene, id, ioshape.radius);
    set_shape_tangents(scene, id, ioshape.tangents);
    set_shape_frames(scene, id, ioshape.instances, ioshape.frame);
    set_shape_emission(scene, id, ioshape.emission, ioshape.emission_tex);
    set_shape_color(scene, id, ioshape.color, ioshape.color_tex);
    set_shape_specular(scene, id, ioshape.specular, ioshape.specular_tex);
    set_shape_ior(scene, id, ioshape.ior);
    set_shape_metallic(scene, id, ioshape.metallic, ioshape.metallic_tex);
    set_shape_transmission(scene, id, ioshape.transmission, ioshape.thin,
        ioshape.trdepth, ioshape.transmission_tex);
    set_shape_roughness(scene, id, ioshape.roughness, ioshape.roughness_tex);
    set_shape_opacity(scene, id, ioshape.opacity, ioshape.opacity_tex);
    set_shape_thin(scene, id, ioshape.thin);
    set_shape_normalmap(scene, id, ioshape.normal_tex);
    set_shape_scattering(scene, id, ioshape.scattering, ioshape.scanisotropy,
        ioshape.scattering_tex);
    ioshape = {};
  }

  for (auto& ioenvironment : ioscene.environments) {
    auto id = add_environment(scene);
    set_environment_frame(scene, id, ioenvironment.frame);
    set_environment_emission(
        scene, id, ioenvironment.emission, ioenvironment.emission_tex);
  }

  ioscene = {};
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
  auto ioscene    = sceneio_model{};
  auto load_timer = print_timed("loading scene");
  load_scene(filename, ioscene);
  print_elapsed(load_timer);

  // add components
  if (validate) {
    auto validate_timer = print_timed("validating");
    auto errors         = scene_validation(ioscene);
    for (auto& error : errors) print_info(error);
    print_elapsed(validate_timer);
  }

  // convert scene
  auto convert_timer = print_timed("converting");
  auto scene         = trace_scene{};
  init_scene(scene, ioscene);
  print_elapsed(convert_timer);

  // build bvh
  auto bvh_timer = print_timed("building bvh");
  init_bvh(scene, params);
  print_elapsed(bvh_timer);

  // init renderer
  auto lights_timer = print_timed("building lights");
  init_lights(scene);
  print_elapsed(lights_timer);

  // fix renderer type if no lights
  if (scene.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, switching to eyelight shader");
    params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  auto state = trace_state{};
  init_state(state, scene, params);
  auto render = image{state.size(), zero4f};

  // render
  for (auto sample = 0; sample < params.samples; sample += batch) {
    auto nsamples    = min(batch, params.samples - sample);
    auto batch_timer = print_timed("rendering samples " +
                                   std::to_string(sample) + "/" +
                                   std::to_string(params.samples));
    render           = trace_samples(state, scene, nsamples, params);
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
