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

#include <map>

// Construct a scene from io
trace_scene make_scene(sceneio_model& ioscene) {
  auto scene = trace_scene{};

  for (auto& iocamera : ioscene.cameras) {
    auto& camera = scene.cameras.emplace_back();
    camera.frame = iocamera.frame;
    camera.film  = iocamera.aspect >= 1
                      ? vec2f{iocamera.film, iocamera.film / iocamera.aspect}
                      : vec2f{iocamera.film * iocamera.aspect, iocamera.film};
    camera.lens     = iocamera.lens;
    camera.focus    = iocamera.focus;
    camera.aperture = iocamera.aperture;
  }

  for (auto& iotexture : ioscene.textures) {
    auto& texture = scene.textures.emplace_back();
    swap(texture.hdr, iotexture.hdr);
    swap(texture.ldr, iotexture.ldr);
  }

  for (auto& iomaterial : ioscene.materials) {
    auto& material            = scene.materials.emplace_back();
    material.emission         = iomaterial.emission;
    material.diffuse          = iomaterial.diffuse;
    material.specular         = iomaterial.specular;
    material.transmission     = iomaterial.transmission;
    material.roughness        = iomaterial.roughness;
    material.opacity          = iomaterial.opacity;
    material.refract          = iomaterial.refract;
    material.volemission      = iomaterial.volemission;
    material.voltransmission  = iomaterial.voltransmission;
    material.volmeanfreepath  = iomaterial.volmeanfreepath;
    material.volscatter       = iomaterial.volscatter;
    material.volscale         = iomaterial.volscale;
    material.volanisotropy    = iomaterial.volanisotropy;
    material.emission_tex     = iomaterial.emission_tex;
    material.diffuse_tex      = iomaterial.diffuse_tex;
    material.specular_tex     = iomaterial.specular_tex;
    material.transmission_tex = iomaterial.transmission_tex;
    material.roughness_tex    = iomaterial.roughness_tex;
    material.opacity_tex      = iomaterial.opacity_tex;
    material.subsurface_tex   = iomaterial.subsurface_tex;
    material.normal_tex       = iomaterial.normal_tex;
  }

  for (auto& ioshape_ : ioscene.shapes) {
    auto tshape = (needs_tesselation(ioscene, ioshape_))
                      ? tesselate_shape(ioscene, ioshape_)
                      : sceneio_shape{};
    auto& ioshape = (needs_tesselation(ioscene, ioshape_)) ? tshape : ioshape_;
    auto& shape   = scene.shapes.emplace_back();
    swap(shape.points, ioshape.points);
    swap(shape.lines, ioshape.lines);
    swap(shape.triangles, ioshape.triangles);
    swap(shape.quads, ioshape.quads);
    swap(shape.quadspos, ioshape.quadspos);
    swap(shape.quadsnorm, ioshape.quadsnorm);
    swap(shape.quadstexcoord, ioshape.quadstexcoord);
    swap(shape.positions, ioshape.positions);
    swap(shape.normals, ioshape.normals);
    swap(shape.texcoords, ioshape.texcoords);
    swap(shape.colors, ioshape.colors);
    swap(shape.radius, ioshape.radius);
    swap(shape.tangents, ioshape.tangents);
    tshape  = {};
    ioshape = {};
  }

  for (auto& ioinstance : ioscene.instances) {
    auto& instance    = scene.instances.emplace_back();
    instance.frame    = ioinstance.frame;
    instance.shape    = ioinstance.shape;
    instance.material = ioinstance.material;
  }

  for (auto& ioenvironment : ioscene.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = ioenvironment.frame;
    environment.emission     = ioenvironment.emission;
    environment.emission_tex = ioenvironment.emission_tex;
  }

  ioscene = {};  // clear
  return scene;
}

int main(int argc, const char* argv[]) {
  // options
  auto params     = trace_params{};
  auto batch      = 16;
  auto save_batch = false;
  auto add_skyenv = false;
  auto validate   = false;
  auto imfilename = "out.hdr"s;
  auto filename   = "scene.json"s;

  // names for enums
  auto sampler_namemap = std::map<string, trace_sampler_type>{};
  for (auto type = 0; type < trace_sampler_names.size(); type++) {
    sampler_namemap[trace_sampler_names[type]] = (trace_sampler_type)type;
  }
  auto falsecolor_namemap = std::map<string, trace_falsecolor_type>{};
  for (auto type = 0; type < trace_falsecolor_names.size(); type++) {
    falsecolor_namemap[trace_falsecolor_names[type]] =
        (trace_falsecolor_type)type;
  }

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
  if (!parse_cli(cli, argc, argv)) exit(1);

  // scene loading
  auto ioscene    = sceneio_model{};
  auto load_timer = print_timed("loading scene");
  if (auto ret = load_scene(filename, ioscene); !ret) {
    print_fatal(ret.error);
  }
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
  auto scene         = make_scene(ioscene);
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
  auto state  = make_state(scene, params);
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
      if (auto ret = save_image(outfilename, render); !ret)
        print_fatal(ret.error);
    }
  }

  // save image
  auto save_timer = print_timed("saving image");
  if (auto ret = save_image(imfilename, render); !ret) print_fatal(ret.error);
  print_elapsed(save_timer);

  // done
  return 0;
}
