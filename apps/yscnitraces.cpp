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
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
#include "ysceneuit.h"
using namespace yocto;

#include <atomic>
#include <future>
#include <map>
#include <thread>

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};  // namespace yocto

// Application scene
struct app_scene {
  // loading options
  string filename  = "scene.yaml";
  string imagename = "out.png";
  string outname   = "out.yaml";
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
  bvh_scene   bvh        = {};
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

// Application state
struct app_state {
  // data
  std::deque<app_scene>    scenes;
  int                      selected = -1;
  std::deque<string>       errors;
  std::deque<app_scene>    loading;
  std::deque<future<void>> load_workers;

  // default options
  load_params    load_prms    = {};
  save_params    save_prms    = {};
  bvh_params     bvh_prms     = {};
  trace_params   trace_prms   = {};
  tonemap_params tonemap_prms = {};
  bool           add_skyenv   = false;
};

void reset_render(app_scene& scene) {
  auto image_size = camera_resolution(
      scene.scene.cameras[scene.trace_prms.camera],
      scene.trace_prms.resolution);
  scene.render.resize(image_size);
  scene.display.resize(image_size);
  scene.render_preview = true;
  scene.render_sample = 0;
  scene.render_region = 0;
  scene.state = make_trace_state(scene.render.size(), scene.trace_prms.seed);
  scene.render_regions = make_regions(
      scene.render.size(), scene.trace_prms.region, true);
}

void load_scene_async(app_state& app, const string& filename) {
  auto& scene        = app.loading.emplace_back();
  scene.filename     = filename;
  scene.imagename    = replace_extension(filename, ".png");
  scene.outname      = replace_extension(filename, ".edited.yaml");
  scene.name         = get_filename(scene.filename);
  scene.load_prms    = app.load_prms;
  scene.save_prms    = app.save_prms;
  scene.trace_prms   = app.trace_prms;
  scene.bvh_prms     = app.bvh_prms;
  scene.tonemap_prms = app.tonemap_prms;
  scene.add_skyenv   = app.add_skyenv;
  app.load_workers.push_back(run_async([&scene]() {
    load_scene(scene.filename, scene.scene);
    make_bvh(scene.bvh, scene.scene, scene.bvh_prms);
    make_trace_lights(scene.lights, scene.scene);
    auto image_size = camera_resolution(
        scene.scene.cameras[scene.trace_prms.camera],
        scene.trace_prms.resolution);
    scene.render.resize(image_size);
    scene.display.resize(image_size);
    scene.name = get_filename(scene.filename) + " [" +
                 std::to_string(scene.render.size().x) + "x" +
                 std::to_string(scene.render.size().y) + " @ 0]";
  }));
}

bool draw_glwidgets_camera(const opengl_window& win, app_scene& scene, int id) {
  auto& scamera = scene.scene.cameras[id];
  auto  camera  = scamera;
  auto  edited  = 0;
  edited += (int)draw_gltextinput(win, "name", camera.name);
  edited += (int)draw_glslider(win, "frame.x", camera.frame.x, -1, 1);
  edited += (int)draw_glslider(win, "frame.y", camera.frame.y, -1, 1);
  edited += (int)draw_glslider(win, "frame.z", camera.frame.z, -1, 1);
  edited += (int)draw_glslider(win, "frame.o", camera.frame.o, -10, 10);
  edited += (int)draw_glcheckbox(win, "ortho", camera.orthographic);
  edited += (int)draw_glslider(win, "lens", camera.lens, 0.01f, 1);
  edited += (int)draw_glslider(win, "film", camera.film, 0.01f, 0.1f);
  edited += (int)draw_glslider(win, "focus", camera.focus, 0.01f, 1000);
  edited += (int)draw_glslider(win, "aperture", camera.aperture, 0, 5);
  auto from         = camera.frame.o,
       to           = camera.frame.o - camera.focus * camera.frame.z;
  auto from_changed = draw_glslider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_glslider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    camera.frame = lookat_frame(from, to, {0, 1, 0});
    camera.focus = length(from - to);
    edited += 1;
  }
  if (edited) {
    // stop_render_async(scene);
    scamera = camera;
    // start_render_async(scene);
  }
  return edited;
}

/// Visit struct elements.
bool draw_glwidgets_texture(
    const opengl_window& win, app_scene& scene, int id) {
  auto& stexture   = scene.scene.textures[id];
  auto  texture    = yocto_texture{};  // copy
  texture.name     = stexture.name;
  texture.filename = stexture.filename;
  auto edited      = 0;
  edited += draw_gltextinput(win, "name", texture.name);
  edited += draw_gltextinput(win, "filename", texture.filename);
  draw_gllabel(
      win, "hdr", "%d x %d", stexture.hdr.size().x, stexture.hdr.size().y);
  draw_gllabel(
      win, "ldr", "%d x %d", stexture.ldr.size().x, stexture.ldr.size().y);
  if (edited) {
    auto reload = texture.filename != stexture.filename;
    // stop_render_async(scene);
    stexture.name     = texture.name;
    stexture.filename = texture.filename;
    if (reload) {
      try {
        if (is_hdr_filename(texture.filename)) {
          load_image(texture.filename, texture.hdr);
        } else {
          load_imageb(texture.filename, texture.ldr);
        }
      } catch (std::exception& e) {
        push_glmessage("cannot load " + texture.filename);
        log_glinfo(win, "cannot load " + texture.filename);
        log_glinfo(win, e.what());
      }
    }
    // start_render_async(scene);
  }
  return edited;
}

bool draw_glwidgets_material(
    const opengl_window& win, app_scene& scene, int id) {
  auto& smaterial = scene.scene.materials[id];
  auto  material  = smaterial;
  auto  edited    = 0;
  edited += draw_gltextinput(win, "name", material.name);
  edited += draw_glhdrcoloredit(win, "emission", material.emission);
  edited += draw_glcoloredit(win, "diffuse", material.diffuse);
  edited += draw_glcoloredit(win, "specular", material.specular);
  edited += draw_glslider(win, "metallic", material.metallic, 0, 1);
  edited += draw_glslider(win, "roughness", material.roughness, 0, 1);
  edited += draw_glcoloredit(win, "coat", material.coat);
  edited += draw_glcoloredit(win, "transmission", material.transmission);
  edited += draw_glcheckbox(win, "refract", material.refract);
  edited += draw_glcoloredit(win, "vol transmission", material.voltransmission);
  edited += draw_glcoloredit(win, "vol meanfreepath", material.volmeanfreepath);
  edited += draw_glcoloredit(win, "vol scatter", material.volscatter);
  edited += draw_glcoloredit(win, "vol emission", material.volemission);
  edited += draw_glslider(win, "vol scale", material.volscale, 0, 1);
  edited += draw_glslider(win, "vol anisotropy", material.volanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", material.opacity, 0, 1);

  edited += draw_glcombobox(
      win, "emission_tex", material.emission_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, scene.scene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, scene.scene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  if (edited) {
    // stop_render_async(scene);
    smaterial = material;
    //start_render_async(scene);
  }
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_scene& scene, int id) {
  auto& sshape   = scene.scene.shapes[id];
  auto  shape    = yocto_shape{};  // copy
  shape.name     = sshape.name;
  shape.filename = sshape.filename;
  auto edited    = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  draw_gllabel(win, "points", "%ld", sshape.points.size());
  draw_gllabel(win, "lines", "%ld", sshape.lines.size());
  draw_gllabel(win, "triangles", "%ld", sshape.triangles.size());
  draw_gllabel(win, "quads", "%ld", sshape.quads.size());
  draw_gllabel(win, "quads pos", "%ld", sshape.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", sshape.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", sshape.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", sshape.positions.size());
  draw_gllabel(win, "norm", "%ld", sshape.normals.size());
  draw_gllabel(win, "texcoord", "%ld", sshape.texcoords.size());
  draw_gllabel(win, "color", "%ld", sshape.colors.size());
  draw_gllabel(win, "radius", "%ld", sshape.radius.size());
  draw_gllabel(win, "tangsp", "%ld", sshape.tangents.size());
  if (edited) {
    // auto reload = sshape.filename != shape.filename;
    // stop_render_async(scene);
    sshape.name     = shape.name;
    sshape.filename = shape.filename;
    // if (reload) {
    //   try {
    //     load_shape(sshape.filename, sshape.points, sshape.lines,
    //         sshape.triangles, sshape.quads, sshape.quadspos, sshape.quadsnorm,
    //         sshape.quadstexcoord, sshape.positions, sshape.normals,
    //         sshape.texcoords, sshape.colors, sshape.radius);
    //   } catch (std::exception& e) {
    //     push_glmessage("cannot load " + shape.filename);
    //     log_glinfo(win, "cannot load " + shape.filename);
    //     log_glinfo(win, e.what());
    //   }
    //   // add bvh refit here
    // }
    // start_render_async(scene);
  }
  return edited;
}

#if 0
inline bool draw_glwidgets_subdiv(const opengl_window& win, app_scene& scene, int id) {
  auto& ssubdiv   = scene.scene.subdiv[id];
  auto  subdiv    = yocto_subdiv{};  // copy
  edited.name             = value.name;
  edited.filename         = value.filename;
  edited.subdivisions     = value.subdivisions;
  edited.catmullclark     = value.catmullclark;
  edited.smooth           = value.smooth;
  edited.facevarying      = value.facevarying;
  edited.shape            = value.shape;
  edited.displacement_tex = value.displacement_tex;
  edited.displacement     = value.displacement;
  auto edited            = 0;
  if (draw_gltextinput(win, "name", edited.name)) updated = true;
  if (draw_gltextinput(win, "filename", edited.filename)) updated = true;
  if (draw_glslider(win, "subdivisions", edited.subdivisions, 0, 10))
    updated = true;
  if (draw_glcheckbox(win, "catmullclark", edited.catmullclark)) updated = true;
  if (draw_glcheckbox(win, "smooth", edited.smooth)) updated = true;
  if (draw_glcheckbox(win, "facevarying", edited.facevarying)) updated = true;
  if (draw_glcombobox(win, "shape", edited.shape, scene.textures, true))
    updated = true;
  if (draw_glcombobox(win, "displacement_tex", edited.displacement_tex,
          scene.textures, true))
    updated = true;
  if (draw_glslider(win, "displacement", edited.displacement, 0, 1))
    updated = true;
  draw_gllabel(win, "points", "%ld", value.points.size());
  draw_gllabel(win, "lines", "%ld", value.lines.size());
  draw_gllabel(win, "triangles", "%ld", value.triangles.size());
  draw_gllabel(win, "quads", "%ld", value.quads.size());
  draw_gllabel(win, "quads pos", "%ld", value.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", value.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", value.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", value.positions.size());
  draw_gllabel(win, "norm", "%ld", value.normals.size());
  draw_gllabel(win, "texcoord", "%ld", value.texcoords.size());
  draw_gllabel(win, "color", "%ld", value.colors.size());
  draw_gllabel(win, "radius", "%ld", value.radius.size());
  if (updated) {
    auto reload = edited.filename != value.filename;
    if (!reload) {
      edited.points        = value.points;
      edited.lines         = value.lines;
      edited.triangles     = value.triangles;
      edited.quads         = value.quads;
      edited.quadspos      = value.quadspos;
      edited.quadsnorm     = value.quadsnorm;
      edited.quadstexcoord = value.quadstexcoord;
      edited.positions     = value.positions;
      edited.normals       = value.normals;
      edited.texcoords     = value.texcoords;
      edited.colors        = value.colors;
      edited.radius        = value.radius;
    }
    edit = {sel.type, sel.index, edited, reload};
  }
  return updated;
}
#endif

#if 0
bool draw_glwidgets_instance(const opengl_window& win, app_scene& scene, int id) {
  auto 
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "name", edited.name)) updated = true;
  if (draw_glslider(win, "frame[0]", edited.frame.x, -1, 1)) updated = true;
  if (draw_glslider(win, "frame[1]", edited.frame.y, -1, 1)) updated = true;
  if (draw_glslider(win, "frame[2]", edited.frame.z, -1, 1)) updated = true;
  if (draw_glslider(win, "frame.o", edited.frame.o, -10, 10)) updated = true;
  if (draw_glcombobox(win, "shape", edited.shape, scene.shapes, true))
    updated = true;
  if (draw_glcombobox(win, "material", edited.material, scene.materials, true))
    updated = true;
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}

bool draw_glwidgets_environment(const opengl_window& win, app_scene& scene, int id) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "name", edited.name)) updated = true;
  if (draw_glslider(win, "frame[0]", edited.frame.x, -1, 1)) updated = true;
  if (draw_glslider(win, "frame[1]", edited.frame.y, -1, 1)) updated = true;
  if (draw_glslider(win, "frame[2]", edited.frame.z, -1, 1)) updated = true;
  if (draw_glslider(win, "frame.o", edited.frame.o, -10, 10)) updated = true;
  if (draw_glhdrcoloredit(win, "emission", edited.emission)) updated = true;
  if (draw_glcombobox(
          win, "emission texture", edited.emission_tex, scene.textures, true))
    updated = true;
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}
#endif 

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app      = *(app_state*)get_gluser_pointer(win);
  auto          scene_ok = !app.scenes.empty() && app.selected >= 0;
  if (!begin_glwidgets_window(win, "yscnitrace")) return;
  draw_glmessages(win);
  if (draw_glfiledialog(
          win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(app, load_path);
  }
  if (draw_glfiledialog(win, "save", save_path, true, get_dirname(save_path),
          get_filename(save_path), "*.yaml;*.obj;*.pbrt")) {
    app.scenes[app.selected].outname = save_path;
    // TODO: support save
    save_path = "";
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    app.scenes[app.selected].imagename = save_path;
    // TODO: support save
    save_path = "";
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save", scene_ok)) {
    save_path = app.scenes[app.selected].outname;
    open_glmodal(win, "save");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save image", scene_ok)) {
    save_path = app.scenes[app.selected].imagename;
    open_glmodal(win, "save image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", scene_ok)) {
    // TODO: support close
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.scenes.empty()) return;
  draw_glcombobox(
      win, "scene", app.selected, (int)app.scenes.size(),
      [&app](int idx) { return app.scenes[idx].name.c_str(); }, false);
  if (scene_ok && begin_glheader(win, "trace")) {
    auto& scn       = app.scenes[app.selected];
    auto  cam_names = vector<string>();
    for (auto& camera : scn.scene.cameras) cam_names.push_back(camera.name);
    auto trace_prms = scn.trace_prms;
    draw_glcombobox(win, "camera", trace_prms.camera, cam_names);
    draw_glslider(win, "resolution", trace_prms.resolution, 180, 4096);
    draw_glslider(win, "nsamples", trace_prms.samples, 16, 4096);
    draw_glcombobox(
        win, "tracer", (int&)trace_prms.sampler, trace_sampler_names);
    draw_glcombobox(win, "false color", (int&)trace_prms.falsecolor,
        trace_falsecolor_names);
    draw_glslider(win, "nbounces", trace_prms.bounces, 1, 128);
    draw_glcheckbox(win, "env hidden", trace_prms.envhidden);
    continue_glline(win);
    draw_glcheckbox(win, "filter", trace_prms.tentfilter);
    draw_glslider(win, "seed", (int&)trace_prms.seed, 0, 1000000);
    draw_glslider(win, "pratio", scn.preview_ratio, 1, 64);
    auto tonemap_prms = scn.tonemap_prms;
    draw_glslider(win, "exposure", tonemap_prms.exposure, -5, 5);
    draw_glcheckbox(win, "filmic", tonemap_prms.filmic);
    continue_glline(win);
    draw_glcheckbox(win, "srgb", tonemap_prms.srgb);
    if (trace_prms != scn.trace_prms) {
      // TODO: support edit
    }
    if (tonemap_prms != scn.tonemap_prms) {
      // TODO: support edit
    }
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "inspect")) {
    auto& scn = app.scenes[app.selected];
    draw_gllabel(win, "scene", get_filename(scn.filename));
    draw_gllabel(win, "filename", scn.filename);
    draw_gllabel(win, "outname", scn.outname);
    draw_gllabel(win, "imagename", scn.imagename);
    draw_gllabel(win, "image", "%d x %d @ %d", scn.render.size().x,
        scn.render.size().y, scn.render_sample);
    draw_glslider(win, "zoom", scn.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", scn.zoom_to_fit);
    continue_glline(win);
    draw_glcheckbox(win, "fps", scn.navigation_fps);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : scn.scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : format_stats(scn.scene)) print_info(stat);
      for (auto stat : format_stats(scn.bvh)) print_info(stat);
    }
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
    draw_gldragger(win, "mouse", ij);
    if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
        ij.y < scn.render.size().y) {
      draw_glcoloredit(win, "pixel", scn.render[{ij.x, ij.y}]);
    } else {
      auto zero4f_ = zero4f;
      draw_glcoloredit(win, "pixel", zero4f_);
    }
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "selection")) {
    auto& scene  = app.scenes[app.selected];
    auto  labels = vector<string>{};
    if (!scene.scene.cameras.empty()) labels.push_back("camera");
    if (!scene.scene.textures.empty()) labels.push_back("texture");
    if (!scene.scene.materials.empty()) labels.push_back("material");
    if (!scene.scene.shapes.empty()) labels.push_back("shape");
    if (!scene.scene.instances.empty()) labels.push_back("instance");
    if (!scene.scene.environments.empty()) labels.push_back("environment");
    if (draw_glcombobox(win, "selection##1", scene.selection.first, labels))
      scene.selection.second = 0;
    if (scene.selection.first == "camera") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.cameras);
      draw_glwidgets_camera(win, scene, scene.selection.second);
    } else if (scene.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.textures);
      draw_glwidgets_texture(win, scene, scene.selection.second);
    } else if (scene.selection.first == "material") {
    } else if (scene.selection.first == "shape") {
    } else if (scene.selection.first == "instance") {
    }
    end_glheader(win);
  }
#if 0
  if (scene_ok && begin_glheader(win, "scene tree")) {
    auto& scn = app.scenes[app.selected];
    draw_glscenetree(win, "", scn.scene, scn.selection, 200);
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "scene object")) {
    auto& scn  = app.scenes[app.selected];
    auto  edit = app_edit{};
    if (draw_glsceneinspector(win, "", scn.scene, scn.selection, edit, 200)) {
      // TODO: support edit
    }
    end_glheader(win);
  }
#endif
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win) {
  auto& app      = *(app_state*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.scenes.empty() && app.selected >= 0) {
    auto& scene = app.scenes.at(app.selected);
    if (!scene.gl_txt || scene.gl_txt.size != scene.display.size())
      init_gltexture(scene.gl_txt, scene.display, false, false, false);
    update_imview(scene.image_center, scene.image_scale, scene.display.size(),
        win_size, scene.zoom_to_fit);
    draw_glimage_background(scene.gl_txt, win_size.x, win_size.y,
        scene.image_center, scene.image_scale);
    set_glblending(true);
    draw_glimage(scene.gl_txt, win_size.x, win_size.y, scene.image_center,
        scene.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void apply_edit(const string& filename, yocto_scene& scene,
    trace_params& trace_prms, tonemap_params& tonemap_prms,
    bool& reload_element, bool& updated_lights, bool& updated_bvh,
    const app_edit& edit) {
  auto& [type, index, data, reload] = edit;

  if (type == typeid(yocto_camera)) {
    scene.cameras[index] = any_cast<yocto_camera>(data);
  } else if (type == typeid(yocto_texture)) {
    scene.textures[index] = any_cast<yocto_texture>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_voltexture)) {
    scene.voltextures[index] = any_cast<yocto_voltexture>(data);
    if (reload) reload_element = true;
  } else if (type == typeid(yocto_shape)) {
    scene.shapes[index] = any_cast<yocto_shape>(data);
    if (reload) {
      reload_element = true;
      updated_bvh    = true;
    }
  } else if (type == typeid(yocto_subdiv)) {
    scene.subdivs[index] = any_cast<yocto_subdiv>(data);
    if (reload) {
      reload_element = true;
      updated_bvh    = true;
    }
  } else if (type == typeid(yocto_material)) {
    auto old_emission      = scene.materials[index].emission;
    scene.materials[index] = any_cast<yocto_material>(data);
    if (old_emission != scene.materials[index].emission) {
      updated_lights = true;
    }
  } else if (type == typeid(yocto_instance)) {
    auto old_instance      = scene.instances[index];
    scene.instances[index] = any_cast<yocto_instance>(data);
    if (old_instance.shape != scene.instances[index].shape ||
        old_instance.frame != scene.instances[index].frame) {
      updated_bvh = true;
    }
  } else if (type == typeid(yocto_environment)) {
    auto old_emission         = scene.materials[index].emission;
    scene.environments[index] = any_cast<yocto_environment>(data);
    if (old_emission != scene.materials[index].emission) {
      updated_lights = true;
    }
  } else if (type == typeid(trace_params)) {
    trace_prms = any_cast<trace_params>(data);
  } else if (type == typeid(tonemap_params)) {
    tonemap_prms = any_cast<tonemap_params>(data);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }
}

// reload an element
void load_element(
    const string& filename, yocto_scene& scene, const app_edit& edit) {
  auto& [type, index, data, reload] = edit;

  if (type == typeid(yocto_texture)) {
    auto& texture = scene.textures[index];
    if (is_hdr_filename(texture.filename)) {
      load_image(get_dirname(filename) + texture.filename, texture.hdr);
    } else {
      load_imageb(get_dirname(filename) + texture.filename, texture.ldr);
    }
  } else if (type == typeid(yocto_voltexture)) {
    auto& texture = scene.voltextures[index];
    load_volume(get_dirname(filename) + texture.filename, texture.vol);
  } else if (type == typeid(yocto_shape)) {
    auto& shape = scene.shapes[index];
    load_shape(get_dirname(filename) + shape.filename, shape.points,
        shape.lines, shape.triangles, shape.quads, shape.quadspos,
        shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
        shape.texcoords, shape.colors, shape.radius, false);
  } else if (type == typeid(yocto_subdiv)) {
    // TODO: this needs more fixing?
    auto& subdiv = scene.subdivs[index];
    load_shape(get_dirname(filename) + subdiv.filename, subdiv.points,
        subdiv.lines, subdiv.triangles, subdiv.quads, subdiv.quadspos,
        subdiv.quadsnorm, subdiv.quadstexcoord, subdiv.positions,
        subdiv.normals, subdiv.texcoords, subdiv.colors, subdiv.radius,
        subdiv.facevarying);
    tesselate_subdiv(scene, scene.subdivs[index]);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }
}

void refit_bvh(const string& filename, yocto_scene& scene, bvh_scene& bvh,
    const bvh_params& bvh_prms, const app_edit& edit) {
  auto& [type, index, data, reload] = edit;

  auto updated_shapes    = vector<int>{};
  auto updated_instances = vector<int>{};
  if (type == typeid(yocto_shape)) {
    updated_shapes.push_back(index);
  } else if (type == typeid(yocto_subdiv)) {
    auto& subdiv = scene.subdivs[index];
    updated_shapes.push_back(subdiv.shape);
  } else if (type == typeid(yocto_instance)) {
    updated_instances.push_back(index);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }

  refit_bvh(bvh, scene, updated_instances, updated_shapes, bvh_prms);
}

void update(const opengl_window& win, app_state& app) {
  for (auto idx = 0; idx < app.load_workers.size(); idx++) {
    if (!is_ready(app.load_workers[idx])) continue;
    try {
      app.load_workers[idx].get();
      app.scenes.push_back(app.loading[idx]);
      reset_render(app.scenes.back());
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load scene " + app.loading[idx].filename);
      log_glinfo(win, "cannot load scene " + app.loading[idx].filename);
      log_glinfo(win, e.what());
      break;
    }
  }
  for (auto& scene : app.scenes) {
    if(scene.render_sample < 0) {
      // rendering preview
      auto preview_prms = scene.trace_prms;
      preview_prms.resolution /= scene.preview_ratio;
      preview_prms.samples = 1;
      auto preview         = trace_image(
          scene.scene, scene.bvh, scene.lights, preview_prms);
      preview = tonemap(preview, scene.tonemap_prms);
      for (auto j = 0; j < scene.display.size().y; j++) {
        for (auto i = 0; i < scene.display.size().x; i++) {
          auto pi = clamp(i / scene.preview_ratio, 0, preview.size().x - 1),
              pj = clamp(j / scene.preview_ratio, 0, preview.size().y - 1);
          scene.display[{i, j}] = preview[{pi, pj}];
        }
      }
      if(!scene.gl_txt || scene.gl_txt.size != scene.display.size()) {
        init_gltexture(scene.gl_txt, scene.display, false, false, false);
      } else {
        update_gltexture(scene.gl_txt, scene.display, false);
      }
    } else if(scene.render_sample < scene.trace_prms.samples) {
      // rendering blocks
      auto num_regions = min(128, scene.render_regions.size() - scene.render_region);
      parallel_for(scene.render_region, scene.render_region + num_regions,
          [&scene](int region_id) {
            trace_region(scene.render, scene.state, scene.scene, scene.bvh,
                scene.lights, scene.render_regions[region_id], 1, scene.trace_prms);
            tonemap(scene.display, scene.render, scene.render_regions[region_id], 
              scene.tonemap_prms);
          });
      scene.render_region += num_regions;
      if(scene.render_region >= scene.render_regions.size()) {
        scene.render_region = 0;
        scene.render_sample += 1;
      }
      update_gltexture(scene.gl_txt, scene.display, false);
    }
  }
#if 0
  // close if needed
  while (!app.scenes.empty()) {
    auto pos = -1;
    for (auto idx = 0; idx < app.scenes.size(); idx++) {
      for (auto& task : app.scenes[idx].task_queue) {
        if (task.type == app_task_type::close_scene) pos = idx;
      }
    }
    if (pos < 0) break;
    app.scenes.erase(app.scenes.begin() + pos);
    app.selected = app.scenes.empty() ? -1 : 0;
  }

  // consume partial results
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (task.type != app_task_type::render_image) continue;
    auto updated = false;
    while (true) {
      auto region = image_region{};
      {
        std::lock_guard guard{task.queuem};
        if (task.queue.empty()) break;
        region = task.queue.front();
        task.queue.pop_front();
      }
      if (region.size() == zero2i) {
        update_gltexture_region(
            scn.gl_txt, scn.preview, {zero2i, scn.preview.size()}, false);
      } else {
        update_gltexture_region(scn.gl_txt, scn.display, region, false);
      }
      updated = true;
    }
    if (updated) {
      scn.render_sample = max(scn.render_sample, (int)task.current);
      scn.name          = get_filename(scn.filename) + "[" +
                 std::to_string(scn.render.size().x) + "x" +
                 std::to_string(scn.render.size().y) + " @ " +
                 std::to_string(scn.render_sample) + "]";
    }
  }

  // remove unneeded tasks
  for (auto& scn : app.scenes) {
    while (scn.task_queue.size() > 1) {
      auto& task = scn.task_queue.at(0);
      auto& next = scn.task_queue.at(1);
      if (task.type == app_task_type::render_image) {
        if (next.type != app_task_type::render_image &&
            next.type != app_task_type::apply_edit)
          break;
        log_glinfo(win, "cancel rendering " + scn.filename);
      } else if (task.type == app_task_type::apply_edit) {
        if (next.type != app_task_type::apply_edit ||
            task.edit.type != next.edit.type ||
            task.edit.index != next.edit.index)
          break;
        log_glinfo(win, "cancel editing " + scn.filename);
      } else {
        break;
      }
      task.stop = true;
      if (task.result.valid()) {
        try {
          task.result.get();
        } catch (...) {
        }
      }
      scn.task_queue.pop_front();
    }
  }

  // apply synchronous edit
  for (auto& scn : app.scenes) {
    while (!scn.task_queue.empty()) {
      auto& task = scn.task_queue.front();
      if (task.type != app_task_type::apply_edit) break;
      log_glinfo(win, "start editing " + scn.filename);
      try {
        scn.render_done     = false;
        auto reload_element = false, update_bvh = false, update_lights = false;
        apply_edit(scn.filename, scn.scene, scn.trace_prms, scn.tonemap_prms,
            reload_element, update_lights, update_bvh, task.edit);
        log_glinfo(win, "done editing " + scn.filename);
        if (reload_element) {
          scn.load_done = false;
          scn.task_queue.emplace_back(app_task_type::load_element, task.edit);
        }
        if (update_bvh) {
          scn.bvh_done = false;
          scn.task_queue.emplace_back(app_task_type::refit_bvh, task.edit);
        }
        if (update_lights) {
          scn.lights_done = false;
          scn.task_queue.emplace_back(app_task_type::init_lights);
        }
        scn.task_queue.emplace_back(app_task_type::render_image);
      } catch (std::exception& e) {
        log_glerror(win, e.what());
        app.errors.push_back("cannot edit " + scn.filename);
      }
      scn.task_queue.pop_front();
    }
  }

  // grab result of finished tasks
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (!task.result.valid()) continue;
    if (task.type == app_task_type::render_image) {
      std::lock_guard guard{task.queuem};
      if (!task.queue.empty()) continue;
    }
    if (task.result.wait_for(std::chrono::nanoseconds(10)) !=
        std::future_status::ready)
      continue;
    switch (task.type) {
      case app_task_type::none: break;
      case app_task_type::close_scene: break;
      case app_task_type::load_scene: {
        try {
          task.result.get();
          scn.load_done  = true;
          scn.image_size = camera_resolution(
              scn.scene.cameras[scn.trace_prms.camera],
              scn.trace_prms.resolution);
          scn.render.resize(scn.image_size);
          scn.display.resize(scn.image_size);
          scn.preview.resize(scn.image_size);
          scn.name = get_filename(scn.filename) + " [" +
                     std::to_string(scn.render.size().x) + "x" +
                     std::to_string(scn.render.size().y) + " @ 0]";
          log_glinfo(win, "done loading " + scn.filename);
          init_gltexture(scn.gl_txt, scn.display, false, false, false);
          scn.task_queue.emplace_back(app_task_type::build_bvh);
          scn.task_queue.emplace_back(app_task_type::init_lights);
          scn.task_queue.emplace_back(app_task_type::render_image);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = get_filename(scn.filename) + " [errpr]";
          app.errors.push_back("cannot load " + scn.filename);
        }
      } break;
      case app_task_type::load_element: {
        try {
          task.result.get();
          scn.load_done = true;
          log_glinfo(win, "done loading element from " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = get_filename(scn.filename) + " [errpr]";
          app.errors.push_back("cannot load element from " + scn.filename);
        }
      } break;
      case app_task_type::build_bvh: {
        try {
          task.result.get();
          scn.bvh_done = true;
          scn.name     = get_filename(scn.filename);
          log_glinfo(win, "done building bvh " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = get_filename(scn.filename) + " [errpr]";
          app.errors.push_back("cannot build bvh " + scn.filename);
        }
      } break;
      case app_task_type::refit_bvh: {
        try {
          task.result.get();
          scn.bvh_done = true;
          scn.name     = get_filename(scn.filename);
          log_glinfo(win, "done refitting bvh " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = get_filename(scn.filename) + " [errpr]";
          app.errors.push_back("cannot refit bvh " + scn.filename);
        }
      } break;
      case app_task_type::init_lights: {
        try {
          task.result.get();
          scn.lights_done = true;
          scn.name        = get_filename(scn.filename);
          log_glinfo(win, "done building lights " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = get_filename(scn.filename) + " [errpr]";
          app.errors.push_back("cannot build lights " + scn.filename);
        }
      } break;
      case app_task_type::save_image: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.imagename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.imagename);
        }
      } break;
      case app_task_type::save_scene: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.outname);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.outname);
        }
      } break;
      case app_task_type::render_image: {
        try {
          task.result.get();
          scn.render_done = true;
          log_glinfo(win, "done rendering " + scn.filename);
          scn.render_sample = scn.trace_prms.samples;
          scn.name          = get_filename(scn.filename) + " [" +
                     std::to_string(scn.render.size().x) + "x" +
                     std::to_string(scn.render.size().y) + " @ " +
                     std::to_string(scn.render_sample) + "]";
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot render " + scn.filename);
        }
      } break;
      case app_task_type::apply_edit: break;
    }
    scn.task_queue.pop_front();
  }
  // schedule tasks not running
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (task.result.valid()) continue;
    task.stop = false;
    switch (task.type) {
      case app_task_type::none: break;
      case app_task_type::close_scene: break;
      case app_task_type::load_scene: {
        log_glinfo(win, "start loading " + scn.filename);
        scn.load_done   = false;
        scn.bvh_done    = false;
        scn.lights_done = false;
        task.result     = std::async(std::launch::async, [&scn]() {
          load_scene(scn.filename, scn.scene, scn.load_prms);
          tesselate_subdivs(scn.scene);
          if (scn.add_skyenv) add_sky(scn.scene);
        });
      } break;
      case app_task_type::load_element: {
        log_glinfo(win, "start loading element for " + scn.filename);
        scn.load_done = false;
        task.result   = std::async(std::launch::async, [&scn, &task]() {
          load_element(scn.filename, scn.scene, task.edit);
        });
      } break;
      case app_task_type::build_bvh: {
        log_glinfo(win, "start building bvh " + scn.filename);
        scn.bvh_done = false;
        task.result  = std::async(std::launch::async,
            [&scn]() { make_bvh(scn.bvh, scn.scene, scn.bvh_prms); });
      } break;
      case app_task_type::refit_bvh: {
        log_glinfo(win, "start refitting bvh " + scn.filename);
        scn.bvh_done = false;
        task.result  = std::async(std::launch::async, [&scn, &task]() {
          refit_bvh(scn.filename, scn.scene, scn.bvh, scn.bvh_prms, task.edit);
        });
      } break;
      case app_task_type::init_lights: {
        log_glinfo(win, "start building lights " + scn.filename);
        scn.lights_done = false;
        task.result     = std::async(std::launch::async,
            [&scn]() { scn.lights = make_trace_lights(scn.scene); });
      } break;
      case app_task_type::save_image: {
        log_glinfo(win, "start saving " + scn.imagename);
        task.result = std::async(std::launch::async, [&scn]() {
          if (is_hdr_filename(scn.imagename)) {
            save_image(scn.imagename, scn.render);
          } else {
            save_imageb(scn.imagename, tonemapb(scn.render, scn.tonemap_prms));
          }
        });
      } break;
      case app_task_type::save_scene: {
        log_glinfo(win, "start saving " + scn.outname);
        task.result = std::async(std::launch::async,
            [&scn]() { save_scene(scn.outname, scn.scene, scn.save_prms); });
      } break;
      case app_task_type::render_image: {
        log_glinfo(win, "start rendering " + scn.filename);
        scn.render_done = false;
        scn.image_size  = camera_resolution(
            scn.scene.cameras[scn.trace_prms.camera],
            scn.trace_prms.resolution);
        if (scn.lights.instances.empty() && scn.lights.environments.empty() &&
            is_sampler_lit(scn.trace_prms)) {
          log_glinfo(win, "no lights presents, switching to eyelight shader");
          scn.trace_prms.sampler = trace_params::sampler_type::eyelight;
        }
        scn.render_sample = 0;
        scn.name          = get_filename(scn.filename) + " [" +
                   std::to_string(scn.render.size().x) + "x" +
                   std::to_string(scn.render.size().y) + " @ " +
                   std::to_string(scn.render_sample) + "]";
        task.result = std::async(std::launch::async, [&scn, &task]() {
          update_app_render(scn.filename, scn.render, scn.display, scn.preview,
              scn.state, scn.scene, scn.lights, scn.bvh, scn.trace_prms,
              scn.tonemap_prms, scn.preview_ratio, &task.stop, task.current,
              task.queue, task.queuem);
        });
        if (scn.render.size() != scn.image_size) {
          scn.render.resize(scn.image_size);
          scn.display.assign(scn.image_size);
          scn.preview.assign(scn.image_size);
          if (scn.gl_txt) {
            delete_gltexture(scn.gl_txt);
            scn.gl_txt = {};
          }
          init_gltexture(scn.gl_txt, scn.display, false, false, false);
        }
      } break;
      case app_task_type::apply_edit: break;
    }
  }
#endif
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  for (auto& path : paths) load_scene_async(app, path);
}

// run ui loop
void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", &app, draw);
  set_drop_glcallback(win, drop_callback);

  // init widgets
  init_glwidgets(win);

  // loop
  auto mouse_pos = zero2f, last_pos = zero2f;
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto alt_down       = get_glalt_key(win);
    auto shift_down     = get_glshift_key(win);
    auto widgets_active = get_glwidgets_active(win);

    // handle mouse and keyboard for navigation
    auto scene_ok = !app.scenes.empty() && app.selected >= 0;
    if (scene_ok && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& scene      = app.scenes[app.selected];
      auto& camera = scene.scene.cameras.at(scene.trace_prms.camera);
      auto  dolly      = 0.0f;
      auto  pan        = zero2f;
      auto  rotate     = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      reset_render(scene);
    }

    // selection
    if (app.selected >= 0 && (mouse_left || mouse_right) && alt_down &&
        !widgets_active) {
      auto& scn = app.scenes[app.selected];
      auto  ij  = get_image_coords(
          mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
      if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
          ij.y < scn.render.size().y) {
        auto& camera = scn.scene.cameras.at(scn.trace_prms.camera);
        auto  ray    = eval_camera(
            camera, ij, scn.render.size(), {0.5f, 0.5f}, zero2f);
        if (auto isec = intersect_bvh(scn.bvh, ray); isec.hit) {
          scn.selection = {"instance", isec.instance};
        }
      }
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
  app.trace_prms.batch = 1;
  auto no_parallel     = false;
  auto filenames       = vector<string>{};

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
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", app.bvh_prms.use_embree,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-flatten/--no-bvh-embree-flatten",
      app.bvh_prms.embree_flatten, "Flatten embree scene");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      app.bvh_prms.embree_compact, "Embree runs in compact memory");
#endif
  add_cli_option(cli, "--add-skyenv", app.add_skyenv, "Add sky envmap");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (no_parallel) {
    app.bvh_prms.noparallel   = true;
    app.load_prms.noparallel  = true;
    app.save_prms.noparallel  = true;
    app.trace_prms.noparallel = true;
  }

  // loading images
  for (auto filename : filenames) load_scene_async(app, filename);

  // run interactive
  run_ui(app);

  // done
  return 0;
}
