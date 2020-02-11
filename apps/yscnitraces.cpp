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
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <future>
#include <memory>
using namespace std;

// Application state
struct app_state {
  // loading options
  string filename  = "scene.yaml";
  string imagename = "out.png";
  string name      = "";

  // options
  trace_params params = {};
  int          pratio = 8;

  // scene
  shared_ptr<trace_scene> scene      = nullptr;
  bool                    add_skyenv = false;

  // rendering state
  shared_ptr<trace_state> state    = nullptr;
  image<vec4f>            render   = {};
  image<vec4f>            display  = {};
  float                   exposure = 0;

  // view scene
  opengl_image        glimage  = {};
  draw_glimage_params glparams = {};

  // editing
  pair<string, int> selection = {"camera", 0};

  // computation
  int          render_sample  = 0;
  atomic<bool> render_stop    = {};
  future<void> render_future  = {};
  int          render_counter = 0;

  ~app_state() {
    render_stop = true;
    if (render_future.valid()) render_future.get();
  }
};

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

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto        futures  = vector<future<void>>{};
  auto        nthreads = thread::hardware_concurrency();
  atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(async(launch::async, [&func, &next_idx, size]() {
      while (true) {
        auto j = next_idx.fetch_add(1);
        if (j >= size.y) break;
        for (auto i = 0; i < size.x; i++) func({i, j});
      }
    }));
  }
  for (auto& f : futures) f.get();
}

void reset_display(shared_ptr<app_state> app) {
  // stop render
  app->render_stop = true;
  if (app->render_future.valid()) app->render_future.get();

  // reset state
  init_state(app->state.get(), app->scene.get(), app->params);
  app->render.resize(app->state->size());
  app->display.resize(app->state->size());

  // render preview
  auto preview_prms = app->params;
  preview_prms.resolution /= app->pratio;
  preview_prms.samples = 1;
  auto preview         = trace_image(app->scene.get(), preview_prms);
  preview              = tonemap_image(preview, app->exposure);
  for (auto j = 0; j < app->display.size().y; j++) {
    for (auto i = 0; i < app->display.size().x; i++) {
      auto pi              = clamp(i / app->pratio, 0, preview.size().x - 1),
           pj              = clamp(j / app->pratio, 0, preview.size().y - 1);
      app->display[{i, j}] = preview[{pi, pj}];
    }
  }

  // start renderer
  app->render_counter = 0;
  app->render_stop    = false;
  app->render_future  = async(launch::async, [app]() {
    for (auto sample = 0; sample < app->params.samples; sample++) {
      if (app->render_stop) return;
      parallel_for(app->render.size(), [app](const vec2i& ij) {
        if (app->render_stop) return;
        app->render[ij] = trace_sample(
            app->state.get(), app->scene.get(), ij, app->params);
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
    }
  });
}

void run_app(int argc, const char* argv[]) {
  // application
  auto app = make_shared<app_state>();

  // parse command line
  auto cli = make_cli("yscnitrace", "progressive path tracing");
  add_cli_option(cli, "--camera", app->params.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app->params.resolution, "Image resolution.");
  add_cli_option(
      cli, "--samples,-s", app->params.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)app->params.sampler, "Tracer type.",
      trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)app->params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", app->params.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", app->params.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", app->params.tentfilter, "Filter image.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", app->params.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(
      cli, "--bvh", (int&)app->params.bvh, "Bvh type", trace_bvh_names);
  add_cli_option(cli, "--add-skyenv", app->add_skyenv, "Add sky envmap");
  add_cli_option(cli, "--output,-o", app->imagename, "Image output", false);
  add_cli_option(cli, "scene", app->filename, "Scene filename", true);
  parse_cli(cli, argc, argv);

  // scene loading
  auto ioscene    = make_shared<sceneio_model>();
  auto load_timer = print_timed("loading scene");
  load_scene(app->filename, ioscene.get());
  print_elapsed(load_timer);

  // conversion
  auto convert_timer = print_timed("converting");
  app->scene         = make_shared<trace_scene>();
  init_scene(app->scene.get(), ioscene.get());
  print_elapsed(convert_timer);

  // cleanup
  ioscene = nullptr;

  // build bvh
  auto bvh_timer = print_timed("building bvh");
  init_bvh(app->scene.get(), app->params);
  print_elapsed(bvh_timer);

  // init renderer
  auto lights_timer = print_timed("building lights");
  init_lights(app->scene.get());
  print_elapsed(lights_timer);

  // fix renderer type if no lights
  if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
    print_info("no lights presents, switching to eyelight shader");
    app->params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  app->state = make_shared<trace_state>();
  init_state(app->state.get(), app->scene.get(), app->params);
  app->render  = image{app->state->size(), zero4f};
  app->display = app->render;
  reset_display(app);

  // window
  auto win = opengl_window{};
  init_glwindow(win, {1280 + 320, 720}, "yscnitraces", false);

  // callbacks
  set_draw_glcallback(
      win, [app](const opengl_window& win, const opengl_input& input) {
        if (!is_initialized(app->glimage)) init_glimage(app->glimage);
        if (!app->render_counter)
          set_glimage(app->glimage, app->display, false, false);
        app->glparams.window      = input.window_size;
        app->glparams.framebuffer = input.framebuffer_viewport;
        update_imview(app->glparams.center, app->glparams.scale,
            app->display.size(), app->glparams.window, app->glparams.fit);
        draw_glimage(app->glimage, app->glparams);
        app->render_counter++;
        if (app->render_counter > 10) app->render_counter = 0;
      });
  set_uiupdate_glcallback(
      win, [app](const opengl_window& win, const opengl_input& input) {
        if ((input.mouse_left || input.mouse_right) && !input.modifier_alt) {
          auto camera = app->scene->cameras.at(app->params.camera);
          auto dolly  = 0.0f;
          auto pan    = zero2f;
          auto rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_right)
            dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) * camera->focus / 200.0f;
          pan.x = -pan.x;
          update_turntable(camera->frame, camera->focus, rotate, dolly, pan);
          reset_display(app);
        }
      });

  // run ui
  run_ui(win);

  // clear
  clear_glwindow(win);
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
