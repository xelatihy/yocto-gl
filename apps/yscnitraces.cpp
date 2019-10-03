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

#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

// Application state
struct app_state {
  // loading options
  string filename  = "app.yaml";
  string imagename = "out.png";
  string name      = "";

  // options
  load_params    load_prms     = {};
  save_params    save_prms     = {};
  bvh_params     bvh_prms      = {};
  trace_params   trace_prms    = {};
  tonemap_params tonemap_prms  = {};
  int            preview_ratio = 8;

  // scene
  yocto_scene scene      = {};
  trace_bvh   bvh        = {};
  bool        add_skyenv = false;

  // rendering state
  trace_lights lights  = {};
  trace_state  state   = {};
  image<vec4f> render  = {};
  image<vec4f> display = {};

  // view scene
  vec2f          image_center   = zero2f;
  float          image_scale    = 1;
  bool           zoom_to_fit    = true;
  bool           navigation_fps = false;
  opengl_texture gl_txt         = {};

  // editing
  pair<string, int> selection = {"camera", 0};

  // computation
  bool                 render_preview = true;
  int                  render_sample  = 0;
  int                  render_region  = 0;
  vector<image_region> render_regions = {};
};

void reset_display(app_state& app) {
  auto image_size = camera_resolution(
      app.scene.cameras[app.trace_prms.camera], app.trace_prms.resolution);
  app.render.resize(image_size);
  app.display.resize(image_size);
  app.render_preview = true;
  app.render_sample  = 0;
  app.render_region  = 0;
  app.state          = make_trace_state(app.render.size(), app.trace_prms.seed);
  app.render_regions = make_image_regions(
      app.render.size(), app.trace_prms.region, true);
}

void draw(const opengl_window& win) {
  auto& app      = *(app_state*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.gl_txt || app.gl_txt.size != app.display.size())
    init_gltexture(app.gl_txt, app.display, false, false, false);
  update_imview(app.image_center, app.image_scale, app.display.size(), win_size,
      app.zoom_to_fit);
  draw_glimage_background(
      app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
  set_glblending(true);
  draw_glimage(
      app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
  set_glblending(false);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_state& app) {
  if (app.render_preview) {
    // rendering preview
    auto preview_prms = app.trace_prms;
    preview_prms.resolution /= app.preview_ratio;
    preview_prms.samples = 1;
    auto preview = trace_image(app.scene, app.bvh, app.lights, preview_prms);
    preview      = tonemap_image(preview, app.tonemap_prms);
    for (auto j = 0; j < app.display.size().y; j++) {
      for (auto i = 0; i < app.display.size().x; i++) {
        auto pi = clamp(i / app.preview_ratio, 0, preview.size().x - 1),
             pj = clamp(j / app.preview_ratio, 0, preview.size().y - 1);
        app.display[{i, j}] = preview[{pi, pj}];
      }
    }
    if (!app.gl_txt || app.gl_txt.size != app.display.size()) {
      init_gltexture(app.gl_txt, app.display, false, false, false);
    } else {
      update_gltexture(app.gl_txt, app.display, false);
    }
    app.render_preview = false;
  } else if (app.render_sample < app.trace_prms.samples) {
    // rendering blocks
    auto num_regions = min(128, app.render_regions.size() - app.render_region);
    parallel_for(app.render_region, app.render_region + num_regions,
        [&app](int region_id) {
          trace_region(app.render, app.state, app.scene, app.bvh, app.lights,
              app.render_regions[region_id], 1, app.trace_prms);
          tonemap_region(app.display, app.render, app.render_regions[region_id],
              app.tonemap_prms);
        });
    if (!app.gl_txt || app.gl_txt.size != app.display.size()) {
      init_gltexture(app.gl_txt, app.display, false, false, false);
    } else {
      for (auto idx = 0; idx < num_regions; idx++)
        update_gltexture_region(app.gl_txt, app.display,
            app.render_regions[app.render_region + idx], false);
    }
    app.render_region += num_regions;
    if (app.render_region >= app.render_regions.size()) {
      app.render_region = 0;
      app.render_sample += 1;
    }
  }
}

// run ui loop
void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", &app, draw);

  // loop
  auto mouse_pos = zero2f, last_pos = zero2f;
  while (!should_glwindow_close(win)) {
    last_pos         = mouse_pos;
    mouse_pos        = get_glmouse_pos(win);
    auto mouse_left  = get_glmouse_left(win);
    auto mouse_right = get_glmouse_right(win);
    auto alt_down    = get_glalt_key(win);
    auto shift_down  = get_glshift_key(win);

    // handle mouse and keyboard for navigation
    if ((mouse_left || mouse_right) && !alt_down) {
      auto& camera = app.scene.cameras.at(app.trace_prms.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      reset_display(app);
    }

    // update
    update(win, app);

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // clear
  delete_glwindow(win);
}

int main(int argc, const char* argv[]) {
  // application
  app_state app{};
  auto      no_parallel = false;

  // parse command line
  auto cli = make_cli("yscnitrace", "progressive path tracing");
  add_cli_option(cli, "--camera", app.trace_prms.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app.trace_prms.resolution, "Image resolution.");
  add_cli_option(
      cli, "--samples,-s", app.trace_prms.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)app.trace_prms.sampler,
      "Tracer type.", trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)app.trace_prms.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", app.trace_prms.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", app.trace_prms.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", app.trace_prms.tentfilter, "Filter image.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", app.trace_prms.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(cli, "--parallel,/--no-parallel", no_parallel,
      "Disable parallel execution.");
  add_cli_option(
      cli, "--exposure,-e", app.tonemap_prms.exposure, "Hdr exposure");
  add_cli_option(
      cli, "--filmic/--no-filmic", app.tonemap_prms.filmic, "Hdr filmic");
  add_cli_option(cli, "--srgb/--no-srgb", app.tonemap_prms.srgb, "Hdr srgb");
  add_cli_option(cli, "--bvh-high-quality/--no-bvh-high-quality",
      app.bvh_prms.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", app.bvh_prms.embree,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      app.bvh_prms.compact, "Embree runs in compact memory");
#endif
  add_cli_option(cli, "--add-skyenv", app.add_skyenv, "Add sky envmap");
  add_cli_option(cli, "--output,-o", app.imagename, "Image output", false);
  add_cli_option(cli, "scene", app.filename, "Scene filename", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (no_parallel) {
    app.bvh_prms.noparallel   = true;
    app.load_prms.noparallel  = true;
    app.save_prms.noparallel  = true;
    app.trace_prms.noparallel = true;
  }

  // scene loading
  try {
    auto timer = print_timed("loading scene");
    load_scene(app.filename, app.scene, app.load_prms);
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // tesselate
  {
    auto timer = print_timed("tesselating");
    tesselate_subdivs(app.scene);
  }

  // add sky
  if (app.add_skyenv) add_sky(app.scene);

  // build bvh
  {
    auto timer = print_timed("building bvh");
    make_bvh(app.bvh, app.scene, app.bvh_prms);
  }

  // init renderer
  {
    auto timer = print_timed("building lights");
    make_trace_lights(app.lights, app.scene);
  }

  // fix renderer type if no lights
  if (app.lights.instances.empty() && app.lights.environments.empty() &&
      is_sampler_lit(app.trace_prms)) {
    print_info("no lights presents, switching to eyelight shader");
    app.trace_prms.sampler = trace_params::sampler_type::eyelight;
  }

  // allocate buffers
  auto image_size = camera_resolution(
      app.scene.cameras[app.trace_prms.camera], app.trace_prms.resolution);
  app.render  = image{image_size, zero4f};
  app.state   = make_trace_state(image_size, app.trace_prms.seed);
  app.display = app.render;
  reset_display(app);

  // run interactive
  run_ui(app);

  // done
  return 0;
}
