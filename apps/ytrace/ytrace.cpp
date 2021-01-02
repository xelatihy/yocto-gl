//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
#include <yocto/yocto_json.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_trace.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

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
    auto camera          = add_camera(scene);
    camera->frame        = iocamera->frame;
    camera->lens         = iocamera->lens;
    camera->aspect       = iocamera->aspect;
    camera->film         = iocamera->film;
    camera->orthographic = iocamera->orthographic;
    camera->aperture     = iocamera->aperture;
    camera->focus        = iocamera->focus;
    camera_map[iocamera] = camera;
  }

  auto texture_map     = unordered_map<sceneio_texture*, trace_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb)
      progress_cb("converting textures", progress.x++, progress.y);
    auto texture           = add_texture(scene);
    texture->hdr           = iotexture->hdr;
    texture->ldr           = iotexture->ldr;
    texture_map[iotexture] = texture;
  }

  auto material_map     = unordered_map<sceneio_material*, trace_material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb)
      progress_cb("converting materials", progress.x++, progress.y);
    auto material              = add_material(scene);
    material->emission         = iomaterial->emission;
    material->color            = iomaterial->color;
    material->specular         = iomaterial->specular;
    material->roughness        = iomaterial->roughness;
    material->metallic         = iomaterial->metallic;
    material->ior              = iomaterial->ior;
    material->spectint         = iomaterial->spectint;
    material->coat             = iomaterial->coat;
    material->transmission     = iomaterial->transmission;
    material->translucency     = iomaterial->translucency;
    material->scattering       = iomaterial->scattering;
    material->scanisotropy     = iomaterial->scanisotropy;
    material->trdepth          = iomaterial->trdepth;
    material->opacity          = iomaterial->opacity;
    material->thin             = iomaterial->thin;
    material->emission_tex     = texture_map.at(iomaterial->emission_tex);
    material->color_tex        = texture_map.at(iomaterial->color_tex);
    material->specular_tex     = texture_map.at(iomaterial->specular_tex);
    material->metallic_tex     = texture_map.at(iomaterial->metallic_tex);
    material->roughness_tex    = texture_map.at(iomaterial->roughness_tex);
    material->transmission_tex = texture_map.at(iomaterial->transmission_tex);
    material->translucency_tex = texture_map.at(iomaterial->translucency_tex);
    material->spectint_tex     = texture_map.at(iomaterial->spectint_tex);
    material->scattering_tex   = texture_map.at(iomaterial->scattering_tex);
    material->coat_tex         = texture_map.at(iomaterial->coat_tex);
    material->opacity_tex      = texture_map.at(iomaterial->opacity_tex);
    material->normal_tex       = texture_map.at(iomaterial->normal_tex);
    material_map[iomaterial]   = material;
  }

  auto shape_map     = unordered_map<sceneio_shape*, trace_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("converting shapes", progress.x++, progress.y);
    auto shape              = add_shape(scene);
    shape->points           = ioshape->points;
    shape->lines            = ioshape->lines;
    shape->triangles        = ioshape->triangles;
    shape->quads            = ioshape->quads;
    shape->quadspos         = ioshape->quadspos;
    shape->quadsnorm        = ioshape->quadsnorm;
    shape->quadstexcoord    = ioshape->quadstexcoord;
    shape->positions        = ioshape->positions;
    shape->normals          = ioshape->normals;
    shape->texcoords        = ioshape->texcoords;
    shape->colors           = ioshape->colors;
    shape->radius           = ioshape->radius;
    shape->tangents         = ioshape->tangents;
    shape->subdivisions     = ioshape->subdivisions;
    shape->catmullclark     = ioshape->catmullclark;
    shape->smooth           = ioshape->smooth;
    shape->displacement     = ioshape->displacement;
    shape->displacement_tex = texture_map.at(ioshape->displacement_tex);
    shape_map[ioshape]      = shape;
  }

  for (auto ioinstance : ioscene->instances) {
    if (progress_cb)
      progress_cb("converting instances", progress.x++, progress.y);
    auto instance      = add_instance(scene);
    instance->frame    = ioinstance->frame;
    instance->shape    = shape_map.at(ioinstance->shape);
    instance->material = material_map.at(ioinstance->material);
  }

  for (auto ioenvironment : ioscene->environments) {
    if (progress_cb)
      progress_cb("converting environments", progress.x++, progress.y);
    auto environment          = add_environment(scene);
    environment->frame        = ioenvironment->frame;
    environment->emission     = ioenvironment->emission;
    environment->emission_tex = texture_map.at(ioenvironment->emission_tex);
  }

  // done
  if (progress_cb) progress_cb("converting done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

// render params
struct render_params : trace_params {
  string scene     = "scene.json";
  string output    = "out.png";
  string camera    = "";
  bool   addsky    = false;
  bool   savebatch = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, render_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.scene, "scene", "Scene filename.", true);
  serialize_property(mode, json, value.output, "output", "Output filename.");
  serialize_property(mode, json, value.camera, "camera", "Camera name.");
  serialize_property(mode, json, value.addsky, "addsky", "Add sky.");
  serialize_property(mode, json, value.savebatch, "savebatch", "Save batch.");
  serialize_value(mode, json, (trace_params&)value, description);
  serialize_clipositionals(mode, json, {"scene"});
  serialize_clialternates(mode, json,
      {{"samples", "s"}, {"bounces", "b"}, {"output", "o"}, {"tracer", "t"}});
}

// convert images
int run_render(const render_params& params) {
  // scene loading
  auto ioscene_guard = std::make_unique<sceneio_scene>();
  auto ioscene       = ioscene_guard.get();
  auto ioerror       = string{};
  if (!load_scene(params.scene, ioscene, ioerror, print_progress))
    return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(ioscene);

  // get camera
  auto iocamera = get_camera(ioscene, params.camera);

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
  auto bvh_guard = std::make_unique<trace_bvh>();
  auto bvh       = bvh_guard.get();
  init_bvh(bvh, scene, params, print_progress);

  // init renderer
  auto lights_guard = std::make_unique<trace_lights>();
  auto lights       = lights_guard.get();
  init_lights(lights, scene, params, print_progress);

  // fix renderer type if no lights
  if (lights->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
  }

  // render
  auto render = trace_image(scene, camera, bvh, lights, params, print_progress,
      [savebatch = params.savebatch, output = params.output](
          const image<vec4f>& render, int sample, int samples) {
        if (!savebatch) return;
        auto ext = "-s" + std::to_string(sample + samples) +
                   path_extension(output);
        auto outfilename = replace_extension(output, ext);
        auto ioerror     = ""s;
        print_progress("save image", sample, samples);
        if (!save_image(outfilename, render, ioerror)) print_fatal(ioerror);
      });

  // save image
  print_progress("save image", 0, 1);
  if (!save_image(params.output, render, ioerror)) return print_fatal(ioerror);
  print_progress("save image", 1, 1);

  // done
  return 0;
}

// convert params
struct view_params : trace_params {
  string scene     = "scene.json";
  string output    = "out.png";
  string camera    = "";
  bool   addsky    = false;
  bool   savebatch = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, view_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.scene, "scene", "Scene filename.", true);
  serialize_property(mode, json, value.output, "output", "Output filename.");
  serialize_property(mode, json, value.camera, "camera", "Camera name.");
  serialize_property(mode, json, value.addsky, "addsky", "Add sky.");
  serialize_property(mode, json, value.savebatch, "savebatch", "Save batch.");
  serialize_value(mode, json, (trace_params&)value, description);
  serialize_clipositionals(mode, json, {"scene"});
  serialize_clialternates(mode, json,
      {{"samples", "s"}, {"bounces", "b"}, {"output", "o"}, {"tracer", "t"}});
}

#ifndef YOCTO_OPENGL

// convert images
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// interactive render
int run_view(const view_params& params) {
  // open viewer
  auto viewer_guard = make_imageviewer("yimage");
  auto viewer       = viewer_guard.get();

  // scene loading
  auto ioscene_guard = std::make_unique<sceneio_scene>();
  auto ioscene       = ioscene_guard.get();
  auto ioerror       = string{};
  if (!load_scene(params.scene, ioscene, ioerror, print_progress))
    return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(ioscene);

  // get camera
  auto iocamera = get_camera(ioscene, params.camera);

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
  auto bvh_guard = std::make_unique<trace_bvh>();
  auto bvh       = bvh_guard.get();
  init_bvh(bvh, scene, params, print_progress);

  // init renderer
  auto lights_guard = std::make_unique<trace_lights>();
  auto lights       = lights_guard.get();
  init_lights(lights, scene, params, print_progress);

  // fix renderer type if no lights
  if (lights->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
  }

  // init state
  auto state_guard = std::make_unique<trace_state>();
  auto state       = state_guard.get();

  // render start
  trace_start(
      state, scene, camera, bvh, lights, params,
      [viewer](const string& message, int sample, int nsamples) {
        set_widget(viewer, "render", "sample", to_json(sample),
            to_schema(sample, "Current sample"));
        print_progress(message, sample, nsamples);
      },
      [viewer](const image<vec4f>& render, int current, int total) {
        set_image(viewer, "render", render);
      });

  // show rendering params
  set_widgets(
      viewer, "render", to_json(params), to_schema(params, "Render params"));

  // set callback
  set_callback(viewer,
      [state, scene, camera, bvh, lights, viewer, &params](const string& name,
          const json_value& uiparams, const gui_input& input) {
        if (name != "render") return;
        if (!uiparams.is_null()) {
          trace_stop(state);
          (view_params&)params = from_json<view_params>(uiparams);
          // show rendering params
          set_widgets(viewer, "render", to_json(params),
              to_schema(params, "Render params"));
          trace_start(
              state, scene, camera, bvh, lights, params,
              [viewer](const string& message, int sample, int nsamples) {
                set_widget(viewer, "render", "sample", to_json(sample),
                    to_schema(sample, "Current sample"));
                print_progress(message, sample, nsamples);
              },
              [viewer](const image<vec4f>& render, int current, int total) {
                set_image(viewer, "render", render);
              });
        } else if ((input.mouse_left || input.mouse_right) &&
                   input.mouse_pos != input.mouse_last) {
          trace_stop(state);
          auto dolly  = 0.0f;
          auto pan    = zero2f;
          auto rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_right)
            dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) * camera->focus / 200.0f;
          pan.x                                  = -pan.x;
          std::tie(camera->frame, camera->focus) = camera_turntable(
              camera->frame, camera->focus, rotate, dolly, pan);
          trace_start(
              state, scene, camera, bvh, lights, params,
              [viewer](const string& message, int sample, int nsamples) {
                set_widget(viewer, "render", "sample", to_json(sample),
                    to_schema(sample, "Current sample"));
                print_progress(message, sample, nsamples);
              },
              [viewer](const image<vec4f>& render, int current, int total) {
                set_image(viewer, "render", render);
              });
        }
      });

  // run view
  run_viewer(viewer);

  // stop
  trace_stop(state);

  // // save image
  // print_progress("save image", 0, 1);
  // if (!save_image(params.output, render, ioerror)) return
  // print_fatal(ioerror); print_progress("save image", 1, 1);

  // done
  return 0;
}

#endif

struct app_params {
  string        command = "render";
  render_params render  = {};
  view_params   view    = {};
};

// Json IO
void serialize_value(json_mode mode, json_value& json, app_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_command(mode, json, value.command, "command", "Command.");
  serialize_property(mode, json, value.render, "render", "Render offline.");
  serialize_property(mode, json, value.view, "view", "Render interactively.");
}

int main(int argc, const char* argv[]) {
  // parse cli
  auto params = app_params{};
  parse_cli(params, "Render images from scenes", argc, argv);

  // dispatch commands
  if (params.command == "render") {
    return run_render(params.render);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else {
    print_fatal("unknown command");
    return 1;
  }
}
