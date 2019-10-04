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
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
using namespace yocto;

#include <map>

int main(int argc, const char* argv[]) {
  // options
  auto load_prms    = load_params{};
  auto bvh_prms     = bvh_params{};
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
  auto sampler_namemap = std::map<string, trace_params::sampler_type>{};
  for (auto type = 0; type < trace_sampler_names.size(); type++) {
    sampler_namemap[trace_sampler_names[type]] =
        (trace_params::sampler_type)type;
  }
  auto falsecolor_namemap = std::map<string, trace_params::falsecolor_type>{};
  for (auto type = 0; type < trace_falsecolor_names.size(); type++) {
    falsecolor_namemap[trace_falsecolor_names[type]] =
        (trace_params::falsecolor_type)type;
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
      bvh_prms.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", bvh_prms.embree,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      bvh_prms.compact, "Embree runs in compact memory");
#endif
  add_cli_option(cli, "--add-skyenv", add_skyenv, "Add sky envmap");
  add_cli_option(cli, "--output-image,-o", imfilename, "Image filename");
  add_cli_option(cli, "--validate", validate, "Validate scene");
  add_cli_option(cli, "--logo/--no-logo", logo, "Whether to append a logo");
  add_cli_option(cli, "scene", filename, "Scene filename", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (noparallel) {
    bvh_prms.noparallel  = true;
    load_prms.noparallel = true;
  }

  // scene loading
  auto scene = yocto_scene{};
  try {
    auto timer = print_timed("loading scene");
    load_scene(filename, scene, load_prms);
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // tesselate
  {
    auto timer = print_timed("tesselating");
    update_tesselation(scene);
  }

  // add components
  if (validate) {
    auto timer = print_timed("validating");
    print_validation(scene);
  }

  // add sky
  if (add_skyenv) add_sky(scene);

  // build bvh
  auto bvh = trace_bvh{};
  {
    auto timer = print_timed("building bvh");
    make_bvh(bvh, scene, bvh_prms);
  }

  // init renderer
  auto lights = trace_lights{};
  {
    auto timer = print_timed("building lights");
    make_trace_lights(lights, scene);
  }

  // fix renderer type if no lights
  if (lights.instances.empty() && lights.environments.empty() &&
      is_sampler_lit(trace_prms)) {
    print_info("no lights presents, switching to eyelight shader");
    trace_prms.sampler = trace_params::sampler_type::eyelight;
  }

  // allocate buffers
  auto image_size = camera_resolution(
      scene.cameras[trace_prms.camera], trace_prms.resolution);
  auto render = image{image_size, zero4f};
  auto state  = make_trace_state(image_size, trace_prms.seed);

  // render
  for (auto sample = 0; sample < trace_prms.samples;
       sample += trace_prms.batch) {
    auto nsamples = min(trace_prms.batch, trace_prms.samples - sample);
    auto timer    = print_timed("rendering samples " + std::to_string(sample) +
                             "/" + std::to_string(trace_prms.samples));
    trace_samples(render, state, scene, bvh, lights, sample, trace_prms);
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
    auto timer = print_timed("saving image");
    if (is_hdr_filename(imfilename)) {
      save_image(imfilename, logo ? add_logo(render) : render);
    } else {
      save_imageb(
          imfilename, logo ? add_logo(tonemap_imageb(render, tonemap_prms))
                           : tonemap_imageb(render, tonemap_prms));
    }
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // done
  return 0;
}
