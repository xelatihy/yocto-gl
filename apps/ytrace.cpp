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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

#include <map>

int main(int argc, char* argv[]) {
  // options
  auto load_prms    = load_params{};
  auto bvh_prms     = bvh_params{};
  auto trace_prms   = trace_params{};
  auto tonemap_prms = tonemap_params{};
  auto all_cameras  = false;
  auto noparallel   = false;
  auto save_batch   = false;
  auto add_skyenv   = false;
  auto validate     = false;
  auto logo         = true;
  auto imfilename   = "out.hdr"s;
  auto filename     = "scene.json"s;

  // names for enums
  auto trace_sampler_type_namemap = std::map<string, trace_sampler_type>{};
  for (auto type = 0; type < trace_sampler_names.size(); type++) {
    trace_sampler_type_namemap[trace_sampler_names[type]] =
        (trace_sampler_type)type;
  }
  auto trace_falsecolor_type_namemap =
      std::map<string, trace_falsecolor_type>{};
  for (auto type = 0; type < trace_falsecolor_names.size(); type++) {
    trace_falsecolor_type_namemap[trace_falsecolor_names[type]] =
        (trace_falsecolor_type)type;
  }

  // parse command line
  auto parser = CLI::App{"Offline path tracing"};
  parser.add_option("--camera", trace_prms.camera, "Camera index.");
  parser.add_flag("--all-cameras", all_cameras, "Render all cameras.");
  parser.add_option(
      "--hres,-R", trace_prms.resolution.x, "Image horizontal resolution.");
  parser.add_option(
      "--vres,-r", trace_prms.resolution.y, "Image vertical resolution.");
  parser.add_option("--samples,-s", trace_prms.samples, "Number of samples.");
  parser.add_option("--tracer,-t", trace_prms.sampler, "Trace type.")
      ->transform(CLI::IsMember(trace_sampler_type_namemap));
  parser
      .add_option(
          "--falsecolor,-F", trace_prms.falsecolor, "Tracer false color type.")
      ->transform(CLI::IsMember(trace_falsecolor_type_namemap));
  parser.add_option(
      "--bounces", trace_prms.bounces, "Maximum number of bounces.");
  parser.add_option("--clamp", trace_prms.clamp, "Final pixel clamping.");
  parser.add_flag("--noparallel", noparallel, "Disable parallel execution.");
  parser.add_option(
      "--seed", trace_prms.seed, "Seed for the random number generators.");
  parser.add_option("--batch,-b", trace_prms.batch, "Samples per batch.");
  parser.add_flag("--env-hidden,!--no-env-hidden", trace_prms.envhidden,
      "Environments are hidden in renderer");
  parser.add_option("--save-batch", save_batch, "Save images progressively");
  parser.add_option("--exposure,-e", tonemap_prms.exposure, "Hdr exposure");
  parser.add_flag("--filmic,!--no-filmic", tonemap_prms.filmic, "Hdr filmic");
  parser.add_flag("--srgb,!--no-srgb", tonemap_prms.srgb, "Hdr srgb");
  parser.add_flag("--bvh-high-quality,!--no-bvh-high-quality",
      bvh_prms.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
  parser.add_flag("--bvh-embree,!--no-bvh-embree", bvh_prms.use_embree,
      "Use Embree ratracer");
  parser.add_flag("--bvh-embree-flatten,!--no-bvh-fembree-latten",
      bvh_prms.embree_flatten, "Flatten BVH scene");
  parser.add_flag("--bvh-embree-compact,!--no-bvh-embree-compact",
      bvh_prms.embree_compact, "Embree runs in compact memory");
#endif
  parser.add_flag("--add-skyenv", add_skyenv, "Add sky envmap");
  parser.add_option("--output-image,-o", imfilename, "Image filename");
  parser.add_flag("--validate", validate, "Validate scene");
  parser.add_flag("--logo,!--no-logo", logo, "Whether to append a logo");
  parser.add_option("scene", filename, "Scene filename", true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }

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
    tesselate_subdivs(scene);
  }

  // add components
  if (validate) {
    auto timer = print_timed("validating");
    print_validation(scene);
  }

  // add sky
  if (add_skyenv) add_sky(scene);

  // build bvh
  auto bvh = bvh_scene{};
  {
    auto timer = print_timed("building bvh");
    build_bvh(bvh, scene, bvh_prms);
  }

  // init renderer
  auto lights = trace_lights{};
  {
    auto timer = print_timed("building lights");
    init_trace_lights(lights, scene);
  }

  // fix renderer type if no lights
  if ((lights.instances.empty() && lights.environments.empty()) &&
      is_sampler_lit(trace_prms)) {
    print_info("no lights presents, switching to eyelight shader");
    trace_prms.sampler = trace_sampler_type::eyelight;
  }

  // cameras to render from
  auto selected_cameras = vector<int>{};
  if (all_cameras) {
    for (auto i = 0; i < scene.cameras.size(); i++) {
      selected_cameras.push_back(i);
    }
  } else {
    selected_cameras.push_back(trace_prms.camera);
  }

  // render all selected cameras
  for (auto camera_id : selected_cameras) {
    // set camera
    trace_prms.camera = camera_id;

    // allocate buffers
    auto image_size = camera_image_size(
        scene.cameras[trace_prms.camera], trace_prms.resolution);
    auto render = image{image_size, zero4f};
    auto state  = trace_state{};
    init_trace_state(state, image_size, trace_prms.seed);

    // render
    for (auto sample = 0; sample < trace_prms.samples;
         sample += trace_prms.batch) {
      auto nsamples = min(trace_prms.batch, trace_prms.samples - sample);
      {
        auto timer = print_timed("rendering cam" +
                                 to_string(trace_prms.camera) + " at " +
                                 pad_left(to_string(sample), 4) + " " +
                                 pad_left(to_string(trace_prms.samples), 4));
        trace_samples(render, state, scene, bvh, lights, sample, trace_prms);
      }
      if (save_batch) {
        auto outfilename = get_noextension(imfilename) + ".cam" +
                           to_string(trace_prms.camera) + ".s" +
                           pad_left(to_string(sample + nsamples), 4, '0') +
                           "." + get_extension(imfilename);
        try {
          if (logo) {
            save_tonemapped_with_logo(outfilename, render, tonemap_prms);
          } else {
            save_tonemapped(outfilename, render, tonemap_prms);
          }
        } catch (const std::exception& e) {
          print_fatal(e.what());
        }
      }
    }

    // save image
    try {
      auto outfilename = imfilename;
      if (all_cameras) {
        outfilename = get_noextension(imfilename) + ".cam" +
                      to_string(trace_prms.camera) + "." +
                      get_extension(imfilename);
      }
      auto timer = print_timed("saving image");
      if (logo) {
        save_tonemapped_with_logo(outfilename, render, tonemap_prms);
      } else {
        save_tonemapped(outfilename, render, tonemap_prms);
      }
    } catch (const std::exception& e) {
      print_fatal(e.what());
    }
  }

  // done
  return 0;
}
