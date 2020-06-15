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
#include <yocto/yocto_trace.h>
using namespace yocto;

#include <map>
#include <memory>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

int main(int argc, const char* argv[]) {
  // options
  auto params         = trace_params{};
  auto save_batch     = false;
  auto add_skyenv     = false;
  auto camera_name    = ""s;
  auto imfilename     = "out.hdr"s;
  auto filename       = "scene.json"s;
  auto feature_images = false;

  // parse command line
  auto cli = make_cli("yscntrace", "Offline path tracing");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(cli, "--resolution,-r", params.resolution, "Image resolution.");
  add_option(cli, "--samples,-s", params.samples, "Number of samples.");
  add_option(
      cli, "--tracer,-t", params.sampler, "Trace type.", trace_sampler_names);
  add_option(cli, "--falsecolor,-F", params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_option(cli, "--bounces,-b", params.bounces, "Maximum number of bounces.");
  add_option(cli, "--clamp", params.clamp, "Final pixel clamping.");
  add_option(cli, "--filter/--no-filter", params.tentfilter, "Filter image.");
  add_option(cli, "--env-hidden/--no-env-hidden", params.envhidden,
      "Environments are hidden in renderer");
  add_option(cli, "--save-batch", save_batch, "Save images progressively");
  add_option(cli, "--bvh", params.bvh, "Bvh type", bvh_names);
  add_option(cli, "--skyenv/--no-skyenv", add_skyenv, "Add sky envmap");
  add_option(cli, "--output-image,-o", imfilename, "Image filename");
  add_option(cli, "scene", filename, "Scene filename", true);
  add_option(cli, "--denoise-features,-d", feature_images,
      "Generate denoise feature images");
  parse_cli(cli, argc, argv);

  // scene loading
  auto scene_guard = std::make_unique<scene_model>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (!load_scene(filename, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // add sky
  if (add_skyenv) add_sky(scene);

  // get camera
  auto camera = get_camera(scene, camera_name);

  // tesselation
  tesselate_shapes(scene, print_progress);

  // build bvh
  init_bvh(scene, params, print_progress);

  // init renderer
  init_lights(scene, print_progress);

  // fix renderer type if no lights
  if (scene->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, switching to eyelight shader");
    params.sampler = trace_sampler_type::eyelight;
  }

  // render
  auto render = trace_image(scene, camera, params, print_progress,
      [save_batch, imfilename](
          const image<vec4f>& render, int sample, int samples) {
        if (!save_batch) return;
        auto ext = "-s" + std::to_string(sample + samples) +
                   sfs::path(imfilename).extension().string();
        auto outfilename =
            sfs::path(imfilename).replace_extension(ext).string();
        auto ioerror = ""s;
        print_progress("save image", sample, samples);
        if (!save_image(outfilename, render, ioerror)) print_fatal(ioerror);
      });

  // save image
  print_progress("save image", 0, 1);
  if (!save_image(imfilename, render, ioerror)) print_fatal(ioerror);
  print_progress("save image", 1, 1);

  if (feature_images) {
    const int   feature_bounces = 5;
    const int   feature_samples = 8;
    std::string feature_ext     = "exr"s;

    auto imext = sfs::path(imfilename).extension();
    if (imext != "hdr" && is_hdr_filename(imext)) feature_ext = imext;

    auto base_name =
        sfs::path(imfilename).filename().replace_extension("").string();

    auto fparams    = params;
    fparams.bounces = feature_bounces;
    fparams.samples = feature_samples;

    // render denoise albedo
    fparams.sampler      = trace_sampler_type::albedo;
    auto albedo          = trace_image(scene, camera, fparams, print_progress);
    auto albedo_filename = sfs::path(imfilename)
                               .replace_filename(base_name + "-albedo")
                               .replace_extension(feature_ext)
                               .string();

    print_progress("save albedo feature", 0, 1);
    if (!save_image(albedo_filename, albedo, ioerror)) print_fatal(ioerror);
    print_progress("save albedo feature", 1, 1);

    // render denoise normals
    fparams.sampler      = trace_sampler_type::normal;
    auto normal          = trace_image(scene, camera, fparams, print_progress);
    auto normal_filename = sfs::path(imfilename)
                               .replace_filename(base_name + "-normal")
                               .replace_extension(feature_ext)
                               .string();

    print_progress("save normal feature", 0, 1);
    if (!save_image(normal_filename, normal, ioerror)) print_fatal(ioerror);
    print_progress("save normal feature", 1, 1);
  }

  // done
  return 0;
}
