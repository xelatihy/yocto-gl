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

#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <future>
#include <memory>
using namespace std;

#include "ext/CLI11.hpp"

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
  shared_ptr<sceneio_model> ioscene    = nullptr;
  shared_ptr<trace_scene>   scene      = nullptr;
  bool                      add_skyenv = false;

  // rendering state
  image<vec4f> render   = {};
  image<vec4f> display  = {};
  float        exposure = 0;

  // view scene
  shared_ptr<opengl_image> glimage  = nullptr;
  draw_glimage_params      glparams = {};

  // computation
  int                           render_sample  = 0;
  int                           render_counter = 0;
  shared_ptr<trace_async_state> render_state   = nullptr;

  ~app_state() { trace_async_stop(render_state); }
};

// construct a scene from io
shared_ptr<trace_scene> make_scene(
    shared_ptr<sceneio_model> ioscene, sceneio_progress print_progress = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->subdivs.size() +
             (int)ioscene->instances.size() + (int)ioscene->objects.size()};

  auto scene = make_trace_scene();

  for (auto iocamera : ioscene->cameras) {
    if (print_progress)
      print_progress("convert camera", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
  }

  auto texture_map =
      unordered_map<shared_ptr<sceneio_texture>, shared_ptr<trace_texture>>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (print_progress)
      print_progress("convert texture", progress.x++, progress.y);
    auto texture = add_texture(scene);
    if (!iotexture->hdr.empty()) {
      set_texture(texture, std::move(iotexture->hdr));
    } else if (!iotexture->ldr.empty()) {
      set_texture(texture, std::move(iotexture->ldr));
    }
    texture_map[iotexture] = texture;
  }

  auto material_map =
      unordered_map<shared_ptr<sceneio_material>, shared_ptr<trace_material>>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (print_progress)
      print_progress("convert material", progress.x++, progress.y);
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

  for (auto iosubdiv : ioscene->subdivs) {
    if (print_progress)
      print_progress("convert subdiv", progress.x++, progress.y);
    tesselate_subdiv(ioscene, iosubdiv);
  }

  auto shape_map =
      unordered_map<shared_ptr<sceneio_shape>, shared_ptr<trace_shape>>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (print_progress)
      print_progress("convert shape", progress.x++, progress.y);
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    set_quads(shape, ioshape->quads);
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_colors(shape, ioshape->colors);
    set_radius(shape, ioshape->radius);
    set_tangents(shape, ioshape->tangents);
    shape_map[ioshape] = shape;
  }

  auto instance_map =
      unordered_map<shared_ptr<sceneio_instance>, shared_ptr<trace_instance>>{};
  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    if (print_progress)
      print_progress("convert instance", progress.x++, progress.y);
    auto instance = add_instance(scene);
    set_frames(instance, ioinstance->frames);
    instance_map[ioinstance] = instance;
  }

  for (auto ioobject : ioscene->objects) {
    if (print_progress)
      print_progress("convert object", progress.x++, progress.y);
    auto object = add_object(scene);
    set_frame(object, ioobject->frame);
    set_shape(object, shape_map.at(ioobject->shape));
    set_material(object, material_map.at(ioobject->material));
    set_instance(object, instance_map.at(ioobject->instance));
  }

  for (auto ioenvironment : ioscene->environments) {
    if (print_progress)
      print_progress("convert environment", progress.x++, progress.y);
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // done
  if (print_progress) print_progress("convert done", progress.x++, progress.y);

  return scene;
}

void reset_display(shared_ptr<app_state> app) {
  // stop render
  trace_async_stop(app->render_state);

  // start render
  app->render_counter = 0;
  app->render_state   = trace_async_start(
      app->scene, app->params, {},
      [app = app.get()](const image<vec4f>& render, int current, int total) {
        if (current > 0) return;
        app->render  = render;
        app->display = tonemap_image(app->render, app->exposure);
      },
      [app = app.get()](
          const image<vec4f>& render, int current, int total, const vec2i& ij) {
        app->render[ij]  = render[ij];
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
}

// progress callback
void print_progress(const string& message, int current, int total) {
  static auto pad = [](const string& str, int n) -> string {
    return string(max(0, n - str.size()), '0') + str;
  };
  static auto pade = [](const string& str, int n) -> string {
    return str + string(max(0, n - str.size()), ' ');
  };
  using clock               = std::chrono::high_resolution_clock;
  static int64_t start_time = 0;
  if (current == 0) start_time = clock::now().time_since_epoch().count();
  auto elapsed = clock::now().time_since_epoch().count() - start_time;
  elapsed /= 1000000;  // millisecs
  auto mins  = pad(to_string(elapsed / 60000), 2);
  auto secs  = pad(to_string((elapsed % 60000) / 1000), 2);
  auto msecs = pad(to_string((elapsed % 60000) % 1000), 3);
  auto n     = (int)(30 * (float)current / (float)total);
  auto bar   = "[" + pade(string(n, '='), 30) + "]";
  auto line  = bar + " " + mins + ":" + secs + "." + msecs + " " +
              pade(message, 30);
  printf("\r%s\r", line.c_str());
  if (current == total) printf("\n");
  fflush(stdout);
}

int run_app(int argc, const char* argv[]) {
  // application
  auto app = make_shared<app_state>();

  // maps for getting param
  auto trace_sampler_map = map<string, trace_sampler_type>{};
  for (auto idx = 0; idx < trace_sampler_names.size(); idx++) {
    trace_sampler_map[trace_sampler_names[idx]] = (trace_sampler_type)idx;
  }
  auto trace_falsecolor_map = map<string, trace_falsecolor_type>{};
  for (auto idx = 0; idx < trace_falsecolor_names.size(); idx++) {
    trace_falsecolor_map[trace_falsecolor_names[idx]] =
        (trace_falsecolor_type)idx;
  }
  auto trace_bvh_map = map<string, trace_bvh_type>{};
  for (auto idx = 0; idx < trace_bvh_names.size(); idx++) {
    trace_bvh_map[trace_bvh_names[idx]] = (trace_bvh_type)idx;
  }

  // parse command line
  auto cli = CLI::App{"progressive path tracing"};
  cli.add_option("--camera", app->params.camera, "Camera index.");
  cli.add_option(
      "--resolution,-r", app->params.resolution, "Image resolution.");
  cli.add_option("--samples,-s", app->params.samples, "Number of samples.");
  cli.add_option("--tracer,-t", app->params.sampler, "Tracer type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_map));
  cli.add_option(
         "--falsecolor,-F", app->params.falsecolor, "Tracer false color type.")
      ->transform(CLI::CheckedTransformer(trace_falsecolor_map));
  cli.add_option(
      "--bounces", app->params.bounces, "Maximum number of bounces.");
  cli.add_option("--clamp", app->params.clamp, "Final pixel clamping.");
  cli.add_flag("--filter", app->params.tentfilter, "Filter image.");
  cli.add_flag("--env-hidden,!--no-env-hidden", app->params.envhidden,
      "Environments are hidden in renderer");
  cli.add_option("--bvh", app->params.bvh, "Bvh type")
      ->transform(CLI::CheckedTransformer(trace_bvh_map));
  cli.add_flag("--add-skyenv", app->add_skyenv, "Add sky envmap");
  cli.add_option("--output,-o", app->imagename, "Image output", false);
  cli.add_option("scene", app->filename, "Scene filename")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // scene loading
  app->ioscene = load_scene(app->filename, print_progress);

  // conversion
  app->scene = make_scene(app->ioscene, print_progress);

  // cleanup
  app->ioscene = nullptr;

  // build bvh
  init_bvh(app->scene, app->params, print_progress);

  // init renderer
  init_lights(app->scene, print_progress);

  // fix renderer type if no lights
  if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
    printf("no lights presents, switching to eyelight shader\n");
    app->params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  reset_display(app);

  // window
  auto win = make_glwindow({1280 + 320, 720}, "yscnitraces", false);

  // callbacks
  set_draw_glcallback(
      win, [app](shared_ptr<opengl_window> win, const opengl_input& input) {
        if (!app->glimage) app->glimage = make_glimage();
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
      win, [app](shared_ptr<opengl_window> win, const opengl_input& input) {
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

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
