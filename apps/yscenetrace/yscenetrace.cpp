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
#include <unordered_map>
using std::unordered_map;

// Construct a scene from io
void init_scene(trace_scene* scene, sceneio_scene* ioscene,
    trace_camera*& camera, sceneio_camera* iocamera,
    progress_callback progress_cb = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->instances.size()};

  auto camera_map     = unordered_map<sceneio_camera*, trace_camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb)
      progress_cb("converting cameras", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film,
        iocamera->orthographic);
    set_focus(camera, iocamera->aperture, iocamera->focus);
    camera_map[iocamera] = camera;
  }

  auto texture_map     = unordered_map<sceneio_texture*, trace_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb)
      progress_cb("converting textures", progress.x++, progress.y);
    auto texture = add_texture(scene);
    if (!iotexture->hdr.empty()) {
      set_texture(texture, iotexture->hdr);
    } else if (!iotexture->ldr.empty()) {
      set_texture(texture, iotexture->ldr);
    }
    texture_map[iotexture] = texture;
  }

  auto material_map     = unordered_map<sceneio_material*, trace_material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb)
      progress_cb("converting materials", progress.x++, progress.y);
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
    set_translucency(material, iomaterial->translucency, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->translucency_tex));
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

  auto shape_map     = unordered_map<sceneio_shape*, trace_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("converting shapes", progress.x++, progress.y);
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    set_quads(shape, ioshape->quads);
    set_fvquads(
        shape, ioshape->quadspos, ioshape->quadsnorm, ioshape->quadstexcoord);
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_colors(shape, ioshape->colors);
    set_radius(shape, ioshape->radius);
    set_tangents(shape, ioshape->tangents);
    set_subdivision(
        shape, ioshape->subdivisions, ioshape->catmullclark, ioshape->smooth);
    shape_map[ioshape] = shape;
  }

  for (auto ioinstance : ioscene->instances) {
    if (progress_cb)
      progress_cb("converting instances", progress.x++, progress.y);
    auto instance = add_instance(scene);
    set_frame(instance, ioinstance->frame);
    set_shape(instance, shape_map.at(ioinstance->shape));
    set_material(instance, material_map.at(ioinstance->material));
  }

  for (auto ioenvironment : ioscene->environments) {
    if (progress_cb)
      progress_cb("converting environments", progress.x++, progress.y);
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // done
  if (progress_cb) progress_cb("converting done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

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
  add_option(cli, "--bvh", params.bvh, "Bvh type", trace_bvh_names);
  add_option(cli, "--skyenv/--no-skyenv", add_skyenv, "Add sky envmap");
  add_option(cli, "--output-image,-o", imfilename, "Image filename");
  add_option(cli, "scene", filename, "Scene filename", true);
  add_option(cli, "--denoise-features,-d", feature_images,
      "Generate denoise feature images");
  parse_cli(cli, argc, argv);

  // scene loading
  auto ioscene_guard = std::make_unique<sceneio_scene>();
  auto ioscene       = ioscene_guard.get();
  auto ioerror       = ""s;
  if (!load_scene(filename, ioscene, ioerror, print_progress))
    print_fatal(ioerror);

  // add sky
  if (add_skyenv) add_sky(ioscene);

  // get camera
  auto iocamera = get_camera(ioscene, camera_name);

  // scene conversion
  auto scene_guard = std::make_unique<trace_scene>();
  auto scene       = scene_guard.get();
  auto camera      = (trace_camera*)nullptr;
  init_scene(scene, ioscene, camera, iocamera);

  // cleanup
  ioscene_guard.reset();

  // tesselation
  tesselate_shapes(scene, print_progress);

  // build bvh
  init_bvh(scene, params, print_progress);

  // init renderer
  auto lights_guard = std::make_unique<trace_lights>();
  auto lights       = lights_guard.get();
  init_lights(lights, scene, print_progress);

  // fix renderer type if no lights
  if (lights->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, switching to eyelight shader");
    params.sampler = trace_sampler_type::eyelight;
  }

  // render
  auto render = trace_image(scene, camera, lights, params, print_progress,
      [save_batch, imfilename](
          const image<vec4f>& render, int sample, int samples) {
        if (!save_batch) return;
        auto ext = "-s" + std::to_string(sample + samples) +
                   path_extension(imfilename);
        auto outfilename = replace_extension(imfilename, ext);
        auto ioerror     = ""s;
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

    auto imext = path_extension(imfilename);
    if (imext != "hdr" && is_hdr_filename(imext)) feature_ext = imext;

    auto fparams    = params;
    fparams.bounces = feature_bounces;
    fparams.samples = feature_samples;

    // render denoise albedo
    fparams.sampler = trace_sampler_type::albedo;
    auto albedo = trace_image(scene, camera, lights, fparams, print_progress);
    auto albedo_filename = replace_extension(
        imfilename, "-albedo" + feature_ext);

    print_progress("save albedo feature", 0, 1);
    if (!save_image(albedo_filename, albedo, ioerror)) print_fatal(ioerror);
    print_progress("save albedo feature", 1, 1);

    // render denoise normals
    fparams.sampler = trace_sampler_type::normal;
    auto normal = trace_image(scene, camera, lights, fparams, print_progress);
    auto normal_filename = replace_extension(
        imfilename, "-normal" + feature_ext);

    print_progress("save normal feature", 0, 1);
    if (!save_image(normal_filename, normal, ioerror)) print_fatal(ioerror);
    print_progress("save normal feature", 1, 1);
  }

  // done
  return 0;
}
