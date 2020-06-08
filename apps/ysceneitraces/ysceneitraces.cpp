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
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>
#include <yocto_gui/yocto_gui.h>
using namespace yocto;

#include <future>
#include <memory>

// Application state
struct app_state {
  // loading options
  string filename  = "scene.yaml";
  string imagename = "out.png";
  string name      = "";

  // options
  trace_params params = {};

  // scene
  scene_model*  scene  = new scene_model{};
  scene_camera* camera = nullptr;

  // rendering state
  image<vec4f> render   = {};
  image<vec4f> display  = {};
  float        exposure = 0;

  // view scene
  ogl_image*       glimage  = new ogl_image{};
  ogl_image_params glparams = {};

  // computation
  int          render_sample  = 0;
  int          render_counter = 0;
  trace_state* render_state   = new trace_state{};

  // status
  std::atomic<int> current = 0;
  std::atomic<int> total   = 0;

  ~app_state() {
    if (render_state) {
      trace_stop(render_state);
      delete render_state;
    }
    if (scene) delete scene;
    if (glimage) delete glimage;
  }
};

void reset_display(app_state* app) {
  // stop render
  trace_stop(app->render_state);

  // start render
  app->render_counter = 0;
  trace_start(
      app->render_state, app->scene, app->camera, app->params,
      [app](const string& message, int sample, int nsamples) {
        app->current = sample;
        app->total   = nsamples;
      },
      [app](const image<vec4f>& render, int current, int total) {
        if (current > 0) return;
        app->render  = render;
        app->display = tonemap_image(app->render, app->exposure);
      },
      [app](
          const image<vec4f>& render, int current, int total, const vec2i& ij) {
        app->render[ij]  = render[ij];
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
}

int main(int argc, const char* argv[]) {
  // application
  auto app_guard = std::make_unique<app_state>();
  auto app       = app_guard.get();

  // command line options
  auto camera_name = ""s;
  auto add_skyenv  = false;

  // parse command line
  auto cli = make_cli("yscnitraces", "progressive path tracing");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(
      cli, "--resolution,-r", app->params.resolution, "Image resolution.");
  add_option(cli, "--samples,-s", app->params.samples, "Number of samples.");
  add_option(cli, "--tracer,-t", app->params.sampler, "Tracer type.",
      trace_sampler_names);
  add_option(cli, "--falsecolor,-F", app->params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_option(
      cli, "--bounces,-b", app->params.bounces, "Maximum number of bounces.");
  add_option(cli, "--clamp", app->params.clamp, "Final pixel clamping.");
  add_option(
      cli, "--filter/--no-filter", app->params.tentfilter, "Filter image.");
  add_option(cli, "--env-hidden/--no-env-hidden", app->params.envhidden,
      "Environments are hidden in renderer");
  add_option(cli, "--bvh", app->params.bvh, "Bvh type", bvh_names);
  add_option(cli, "--skyenv/--no-skyenv", add_skyenv, "Add sky envmap");
  add_option(cli, "--output,-o", app->imagename, "Image output");
  add_option(cli, "scene", app->filename, "Scene filename", true);
  parse_cli(cli, argc, argv);

  // scene loading
  auto ioerror = ""s;
  if (!load_scene(app->filename, app->scene, ioerror, print_progress))
    print_fatal(ioerror);

  // add sky
  if (add_skyenv) add_sky(app->scene);

  // get camera
  app->camera = get_camera(app->scene, camera_name);

  // tesselation
  tesselate_shapes(app->scene, print_progress);

  // build bvh
  init_bvh(app->scene, app->params, print_progress);

  // init renderer
  init_lights(app->scene, print_progress);

  // fix renderer type if no lights
  if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
    print_info("no lights presents, switching to eyelight shader");
    app->params.sampler = trace_sampler_type::eyelight;
  }

  // allocate buffers
  reset_display(app);

  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [app](gui_window* win, const gui_input& input) {
    clear_image(app->glimage);
  };
  callbacks.draw_cb = [app](gui_window* win, const gui_input& input) {
    if (!is_initialized(app->glimage)) init_image(app->glimage);
    if (!app->render_counter)
      set_image(app->glimage, app->display, false, false);
    app->glparams.window      = input.window_size;
    app->glparams.framebuffer = input.framebuffer_viewport;
    update_imview(app->glparams.center, app->glparams.scale,
        app->display.size(), app->glparams.window, app->glparams.fit);
    draw_image(app->glimage, app->glparams);
    app->render_counter++;
    if (app->render_counter > 10) app->render_counter = 0;
  };
  callbacks.widgets_cb = [app](gui_window* win, const gui_input& input) {
    auto  edited  = 0;
    auto& tparams = app->params;
    draw_progressbar(win, "render", app->current, app->total);
    edited += draw_combobox(win, "camera", app->camera, app->scene->cameras);
    edited += draw_slider(win, "resolution", tparams.resolution, 180, 4096);
    edited += draw_slider(win, "nsamples", tparams.samples, 16, 4096);
    edited += draw_combobox(
        win, "tracer", (int&)tparams.sampler, trace_sampler_names);
    edited += draw_combobox(
        win, "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
    edited += draw_slider(win, "nbounces", tparams.bounces, 1, 128);
    edited += draw_checkbox(win, "envhidden", tparams.envhidden);
    continue_line(win);
    edited += draw_checkbox(win, "filter", tparams.tentfilter);
    edited += draw_slider(win, "seed", (int&)tparams.seed, 0, 1000000);
    edited += draw_slider(win, "pratio", tparams.pratio, 1, 64);
    edited += draw_slider(win, "exposure", app->exposure, -5, 5);
    if (edited) reset_display(app);
  };
  callbacks.char_cb = [app](gui_window* win, unsigned int key,
                          const gui_input& input) {
    switch (key) {
      case 'c': {
        auto ncameras = (int)app->scene->cameras.size();
        for (auto pos = 0; pos < ncameras; pos++) {
          if (app->scene->cameras[pos] == app->camera) {
            app->camera = app->scene->cameras[(pos + 1) % ncameras];
            reset_display(app);
            break;
          }
        }
      } break;
      case 'f':
        app->params.sampler = trace_sampler_type::falsecolor;
        reset_display(app);
        break;
      case 'p':
        app->params.sampler = trace_sampler_type::path;
        reset_display(app);
        break;
      case 'F':
        app->params.falsecolor = (trace_falsecolor_type)(
            ((int)app->params.falsecolor + 1) %
            (int)trace_sampler_names.size());
        reset_display(app);
        break;
    }
  };
  callbacks.uiupdate_cb = [app](gui_window* win, const gui_input& input) {
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * app->camera->focus /
              200.0f;
      pan.x = -pan.x;
      update_turntable(
          app->camera->frame, app->camera->focus, rotate, dolly, pan);
      reset_display(app);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yscnitraces", callbacks);

  // done
  return 0;
}
