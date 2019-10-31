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
trace_scene make_scene(const scene_model& ioscene) {
  auto tesselate = [](const scene_model& ioscene, const scene_shape& shape) {
    if (!shape.subdivisions && !shape.displacement) return scene_shape{};
    auto subdiv = shape;
    if (subdiv.subdivisions) subdiv = subdivide_shape(subdiv);
    if (subdiv.displacement && subdiv.displacement_tex >= 0)
      subdiv = displace_shape(ioscene, subdiv);
    return subdiv;
  };

  auto scene = trace_scene{};

  for (auto& iocamera : ioscene.cameras) {
    auto& camera = scene.cameras.emplace_back();
    camera.frame = iocamera.frame;
    camera.film  = iocamera.aspect >= 1
                      ? vec2f{iocamera.film, iocamera.film / iocamera.aspect}
                      : vec2f{iocamera.film / iocamera.aspect, iocamera.film};
    camera.lens     = iocamera.lens;
    camera.focus    = iocamera.focus;
    camera.aperture = iocamera.aperture;
  }

  for (auto& iotexture : ioscene.textures) {
    auto& texture = scene.textures.emplace_back();
    texture.hdr   = iotexture.hdr;
    texture.ldr   = iotexture.ldr;
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
  }

  for (auto& ioshape_ : ioscene.shapes) {
    auto& ioshape = (ioshape_.subdivisions || ioshape_.displacement)
                        ? tesselate(ioscene, ioshape_)
                        : ioshape_;
    auto& shape         = scene.shapes.emplace_back();
    shape.points        = ioshape.points;
    shape.lines         = ioshape.lines;
    shape.triangles     = ioshape.triangles;
    shape.quads         = ioshape.quads;
    shape.quadspos      = ioshape.quadspos;
    shape.quadsnorm     = ioshape.quadsnorm;
    shape.quadstexcoord = ioshape.quadstexcoord;
    shape.positions     = ioshape.positions;
    shape.normals       = ioshape.normals;
    shape.texcoords     = ioshape.texcoords;
    shape.colors        = ioshape.colors;
    shape.radius        = ioshape.radius;
    shape.tangents      = ioshape.tangents;
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

  return scene;
}

int main(int argc, const char* argv[]) {
  // options
  auto load_prms    = load_params{};
  auto trace_prms   = trace_params{};
  auto tonemap_prms = tonemap_params{};
  auto noparallel   = false;
  auto save_batch   = false;
  auto add_skyenv   = false;
  auto validate     = false;
  auto logo         = false;
  auto imfilename   = "out.hdr"s;
  auto filename     = "scene.json"s;

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
  add_cli_option(cli, "--camera", trace_prms.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", trace_prms.resolution, "Image resolution.");
  add_cli_option(cli, "--samples,-s", trace_prms.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)trace_prms.sampler, "Trace type.",
      trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)trace_prms.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", trace_prms.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", trace_prms.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", trace_prms.tentfilter, "Filter image.");
  add_cli_option(
      cli, "--noparallel", noparallel, "Disable parallel execution.");
  add_cli_option(cli, "--batch,-b", trace_prms.batch, "Samples per batch.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", trace_prms.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(cli, "--save-batch", save_batch, "Save images progressively");
  add_cli_option(cli, "--exposure,-e", tonemap_prms.exposure, "Hdr exposure");
  add_cli_option(
      cli, "--filmic/--no-filmic", tonemap_prms.filmic, "Hdr filmic");
  add_cli_option(cli, "--srgb/--no-srgb", tonemap_prms.srgb, "Hdr srgb");
  add_cli_option(cli, "--bvh-high-quality/--no-bvh-high-quality",
      trace_prms.highquality_bvh, "Use high quality bvh mode");
#if YOCTO_EMBREE
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", trace_prms.embree_bvh,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      trace_prms.compact_bvh, "Embree runs in compact memory");
#endif
  add_cli_option(cli, "--add-skyenv", add_skyenv, "Add sky envmap");
  add_cli_option(cli, "--output-image,-o", imfilename, "Image filename");
  add_cli_option(cli, "--validate", validate, "Validate scene");
  add_cli_option(cli, "--logo/--no-logo", logo, "Whether to append a logo");
  add_cli_option(cli, "scene", filename, "Scene filename", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (noparallel) {
    load_prms.noparallel = true;
  }

  // scene loading
  auto ioscene    = scene_model{};
  auto load_timer = print_timed("loading scene");
  if (auto ret = load_scene(filename, ioscene, load_prms); !ret) {
    print_fatal(ret.error);
  }
  print_elapsed(load_timer);

  // add components
  if (validate) {
    auto validate_timer = print_timed("validating");
    auto errors         = format_validation(ioscene);
    for (auto& error : errors) print_info(error);
    print_elapsed(validate_timer);
  }

  // convert scene
  auto convert_timer = print_timed("converting");
  auto scene         = make_scene(ioscene);
  print_elapsed(convert_timer);

  // build bvh
  auto bvh_timer = print_timed("building bvh");
  init_bvh(scene, trace_prms);
  print_elapsed(bvh_timer);

  // init renderer
  auto lights_timer = print_timed("building lights");
  init_lights(scene);
  print_elapsed(lights_timer);

  // fix renderer type if no lights
  if (scene.lights.empty() && is_sampler_lit(trace_prms)) {
    print_info("no lights presents, switching to eyelight shader");
    trace_prms.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  auto state  = make_state(scene, trace_prms);
  auto render = image{state.size(), zero4f};

  // render
  for (auto sample = 0; sample < trace_prms.samples;
       sample += trace_prms.batch) {
    auto nsamples    = min(trace_prms.batch, trace_prms.samples - sample);
    auto batch_timer = print_timed("rendering samples " +
                                   std::to_string(sample) + "/" +
                                   std::to_string(trace_prms.samples));
    trace_samples(render, state, scene, trace_prms);
    print_elapsed(batch_timer);
    if (save_batch) {
      auto outfilename = replace_extension(imfilename,
          "-s" + std::to_string(sample + nsamples) + get_extension(imfilename));
      try {
        if (is_hdr_filename(outfilename)) {
          save_image(outfilename, logo ? add_logo(render) : render);
        } else {
          save_imageb(
              outfilename, logo ? add_logo(tonemap_imageb(render, tonemap_prms))
                                : tonemap_imageb(render, tonemap_prms));
        }
      } catch (const std::exception& e) {
        print_fatal(e.what());
      }
    }
  }

  // save image
  try {
    auto save_timer = print_timed("saving image");
    if (is_hdr_filename(imfilename)) {
      save_image(imfilename, logo ? add_logo(render) : render);
    } else {
      save_imageb(
          imfilename, logo ? add_logo(tonemap_imageb(render, tonemap_prms))
                           : tonemap_imageb(render, tonemap_prms));
    }
    print_elapsed(save_timer);
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // done
  return 0;
}
