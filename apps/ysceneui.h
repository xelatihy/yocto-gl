//
// Utilities to display a scene graph using ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_SCENEUI_H_
#define _YOCTO_SCENEUI_H_

#include "../yocto/yocto_scene.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <any>
#include <typeindex>
#include <unordered_map>
using std::any;
using std::any_cast;
using std::type_index;

struct app_selection {
  type_index type  = typeid(void);
  int        index = -1;
};

struct app_edit {
  type_index type   = typeid(void);
  int        index  = -1;
  any        data   = {};
  bool       reload = false;
};

inline const unordered_map<int, string>& animation_type_names() {
  static auto names = unordered_map<int, string>{
      {(int)yocto_animation::type_t::linear, "linear"},
      {(int)yocto_animation::type_t::step, "step"},
      {(int)yocto_animation::type_t::bezier, "bezier"},
  };
  return names;
}

template <typename T>
inline void draw_glscenetree(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T>& vals, app_selection& sel);

template <typename T>
inline void draw_glscenetree(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T*>& vals, app_selection& sel);

inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_camera& value, app_selection& sel) {}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_texture& value, app_selection& sel) {}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_voltexture& value, app_selection& sel) {}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_material& value, app_selection& sel) {
  draw_glscenetree(
      win, "emission", scene, value.emission_tex, scene.textures, sel);
  draw_glscenetree(
      win, "diffuse", scene, value.diffuse_tex, scene.textures, sel);
  draw_glscenetree(
      win, "metallic", scene, value.metallic_tex, scene.textures, sel);
  draw_glscenetree(
      win, "specular", scene, value.specular_tex, scene.textures, sel);
  draw_glscenetree(win, "normal", scene, value.normal_tex, scene.textures, sel);
}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_shape& value, app_selection& sel) {}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_subdiv& value, app_selection& sel) {
  draw_glscenetree(win, "shapes", scene, value.shape, scene.shapes, sel);
  draw_glscenetree(
      win, "displament", scene, value.displacement_tex, scene.textures, sel);
}

inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_instance& value, app_selection& sel) {
  draw_glscenetree(win, "shape", scene, value.shape, scene.shapes, sel);
  draw_glscenetree(
      win, "material", scene, value.material, scene.materials, sel);
}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_environment& value, app_selection& sel) {
  draw_glscenetree(
      win, "emission", scene, value.emission_tex, scene.textures, sel);
}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_scene_node& value, app_selection& sel) {
  draw_glscenetree(
      win, "instance", scene, value.instance, scene.instances, sel);
  draw_glscenetree(win, "camera", scene, value.camera, scene.cameras, sel);
  draw_glscenetree(
      win, "environment", scene, value.environment, scene.environments, sel);
  draw_glscenetree(win, "parent", scene, value.parent, scene.nodes, sel);
  auto cid = 0;
  for (auto ch : value.children) {
    draw_glscenetree(
        win, "child" + std::to_string(cid++), scene, ch, scene.nodes, sel);
  }
}
inline void draw_glscenetree_rec(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, const yocto_animation& value, app_selection& sel) {
  auto tid = 0;
  for (auto tg : value.targets) {
    draw_glscenetree(
        win, "target" + std::to_string(tid++), scene, tg, scene.nodes, sel);
  }
}

template <typename T>
inline void draw_glscenetree(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T>& vals, app_selection& sel) {
  if (index < 0) return;
  auto lbl = vals[index].uri;
  if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index].uri;
  auto selected = sel.type == type_index(typeid(T)) && sel.index == index;
  if (begin_glselectabletreenode(win, lbl.c_str(), selected)) {
    draw_glscenetree_rec(win, lbl_, scene, vals[index], sel);
    end_gltreenode(win);
  }
  if (selected) sel = {type_index(typeid(T)), index};
}

template <typename T>
inline void draw_glscenetree(const opengl_window& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T*>& vals, app_selection& sel) {
  if (index < 0) return;
  auto lbl = vals[index]->name;
  if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index]->name;
  auto selected = sel.type == type_index(typeid(T)) && sel.index == index;
  if (begin_glselectabletreenode(win, lbl.c_str(), selected)) {
    draw_glscenetree_rec(win, lbl_, scene, vals[index], sel);
    end_gltreenode(win);
  }
  if (selected) sel = {type_index(typeid(T)), index};
}

inline void draw_glscenetree(
    const opengl_window& win, yocto_scene& scene, app_selection& sel) {
  if (!empty(scene.cameras) && begin_gltreenode(win, "cameras")) {
    for (auto v = 0; v < scene.cameras.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.cameras, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.shapes) && begin_gltreenode(win, "shapes")) {
    for (auto v = 0; v < scene.shapes.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.shapes, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.instances) && begin_gltreenode(win, "instances")) {
    for (auto v = 0; v < scene.instances.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.instances, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.materials) && begin_gltreenode(win, "materials")) {
    for (auto v = 0; v < scene.materials.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.materials, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.textures) && begin_gltreenode(win, "textures")) {
    for (auto v = 0; v < scene.textures.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.textures, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.environments) && begin_gltreenode(win, "environments")) {
    for (auto v = 0; v < scene.environments.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.environments, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.nodes) && begin_gltreenode(win, "nodes")) {
    for (auto v = 0; v < scene.nodes.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.nodes, sel);
    end_gltreenode(win);
  }
  if (!empty(scene.animations) && begin_gltreenode(win, "animations")) {
    for (auto v = 0; v < scene.animations.size(); v++)
      draw_glscenetree(win, "", scene, v, scene.animations, sel);
    end_gltreenode(win);
  }
}

/// Visit struct elements.
inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_camera& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  if (draw_glslider(win, "frame.x", edited.frame.x, -1, 1)) updated = true;
  if (draw_glslider(win, "frame.y", edited.frame.y, -1, 1)) updated = true;
  if (draw_glslider(win, "frame.z", edited.frame.z, -1, 1)) updated = true;
  if (draw_glslider(win, "frame.o", edited.frame.o, -10, 10)) updated = true;
  if (draw_glcheckbox(win, "ortho", edited.orthographic)) updated = true;
  if (draw_glslider(win, "lens", edited.lens, 0.01f, 1)) updated = true;
  if (draw_glslider(win, "film", edited.film, 0.01f, 0.1f)) updated = true;
  if (draw_glslider(win, "focus", edited.focus, 0.01f, 1000)) updated = true;
  if (draw_glslider(win, "aperture", edited.aperture, 0, 5)) updated = true;
  auto from = edited.frame.o,
       to   = edited.frame.o - edited.focus * edited.frame.z;
  draw_glslider(win, "!!from", from, -10, 10);
  draw_glslider(win, "!!to", to, -10, 10);
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}

/// Visit struct elements.
inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_texture& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = yocto_texture{};
  edited.uri   = value.uri;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  draw_gllabel(win, "hdr", "%d x %d", value.hdr.size().x, value.hdr.size().y);
  draw_gllabel(win, "ldr", "%d x %d", value.ldr.size().x, value.ldr.size().y);
  if (updated) {
    auto reload = edited.uri != value.uri;
    if (!reload) {
      edited.hdr = value.hdr;
      edited.ldr = value.ldr;
    }
    edit = {sel.type, sel.index, edited, reload};
  }
  return updated;
}

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_voltexture& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = yocto_voltexture{};
  edited.uri   = value.uri;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  draw_gllabel(win, "voxel_data", "%d x %d x %d", value.vol.size().x,
      value.vol.size().y, value.vol.size().z);
  if (updated) {
    auto reload = edited.uri != value.uri;
    if (!reload) {
      edited.vol = value.vol;
    }
    edit = {sel.type, sel.index, edited, reload};
  }
  return updated;
}

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_material& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  if (draw_glhdrcoloredit(win, "emission", edited.emission)) updated = true;
  if (draw_glcoloredit(win, "diffuse", edited.diffuse)) updated = true;
  if (draw_glcoloredit(win, "specular", edited.specular)) updated = true;
  if (draw_glslider(win, "metallic", edited.metallic, 0, 1)) updated = true;
  if (draw_glslider(win, "roughness", edited.roughness, 0, 1)) updated = true;
  if (draw_glcoloredit(win, "coat", edited.coat)) updated = true;
  if (draw_glcoloredit(win, "transmission", edited.transmission))
    updated = true;
  if (draw_glcoloredit(win, "vol transmission", edited.voltransmission))
    updated = true;
  if (draw_glcoloredit(win, "vol scatter", edited.volscatter)) updated = true;
  if (draw_glcoloredit(win, "vol emission", edited.volemission)) updated = true;
  if (draw_glslider(win, "vol scale", edited.volscale, 0, 1)) updated = true;
  if (draw_glslider(win, "vol anisotropy", edited.volanisotropy, -1, 1))
    updated = true;
  if (draw_glslider(win, "opacity", edited.opacity, 0, 1)) updated = true;
  if (draw_glcheckbox(win, "thin", edited.thin)) updated = true;

  if (draw_glcombobox(
          win, "emission_tex", edited.emission_tex, scene.textures, true))
    updated = true;
  if (draw_glcombobox(
          win, "diffuse_tex", edited.diffuse_tex, scene.textures, true))
    updated = true;
  if (draw_glcombobox(
          win, "metallic_tex", edited.metallic_tex, scene.textures, true))
    updated = true;
  if (draw_glcombobox(
          win, "specular_tex", edited.specular_tex, scene.textures, true))
    updated = true;
  if (draw_glcombobox(win, "transmission_tex", edited.transmission_tex,
          scene.textures, true))
    updated = true;
  if (draw_glcombobox(
          win, "roughness_tex", edited.roughness_tex, scene.textures, true))
    updated = true;
  if (draw_glcombobox(
          win, "normal_tex", edited.normal_tex, scene.textures, true))
    updated = true;
  if (draw_glcheckbox(win, "glTF textures", edited.gltf_textures))
    updated = true;
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_shape& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = yocto_shape{};
  edited.uri   = value.uri;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
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
  draw_gllabel(win, "tangsp", "%ld", value.tangents.size());
  if (updated) {
    auto reload = edited.uri != value.uri;
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
      edited.tangents      = value.tangents;
    }
    edit = {sel.type, sel.index, edited, reload};
  }
  return updated;
  return updated;
}

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_subdiv& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited             = yocto_subdiv{};
  edited.uri              = value.uri;
  edited.subdivisions     = value.subdivisions;
  edited.catmullclark     = value.catmullclark;
  edited.smooth           = value.smooth;
  edited.facevarying      = value.facevarying;
  edited.shape            = value.shape;
  edited.displacement_tex = value.displacement_tex;
  edited.displacement     = value.displacement;
  auto updated            = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
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
    auto reload = edited.uri != value.uri;
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

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_instance& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
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

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_environment& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
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

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_scene_node& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  if (draw_glcombobox(win, "parent", edited.parent, scene.nodes, true))
    updated = true;
  if (draw_glslider(win, "local[0]", edited.local.x, -1, 1)) updated = true;
  if (draw_glslider(win, "local[1]", edited.local.y, -1, 1)) updated = true;
  if (draw_glslider(win, "local[2]", edited.local.z, -1, 1)) updated = true;
  if (draw_glslider(win, "local.o", edited.local.o, -10, 10)) updated = true;
  if (draw_glslider(win, "translation", edited.translation, -10, 10))
    updated = true;
  if (draw_glslider(win, "rotation", edited.rotation, -1, 1)) updated = true;
  if (draw_glslider(win, "scale", edited.scale, 0, 10)) updated = true;
  if (draw_glcombobox(win, "camera", edited.camera, scene.cameras, true))
    updated = true;
  if (draw_glcombobox(win, "instance", edited.instance, scene.instances, true))
    updated = true;
  if (draw_glcombobox(
          win, "environment", edited.environment, scene.environments, true))
    updated = true;
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}

inline bool draw_glsceneinspector(const opengl_window& win,
    const yocto_animation& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
  auto edited  = value;
  auto updated = false;
  if (draw_gltextinput(win, "uri", edited.uri)) updated = true;
  if (draw_gltextinput(win, "path", edited.filename)) updated = true;
  if (draw_gltextinput(win, "group", edited.group)) updated = true;
  // if(draw_glcombobox(win, "type", &value.type,
  // animation_type_names())) updated = false;
  draw_gllabel(win, "times", "%ld", edited.times.size());
  draw_gllabel(win, "translation", "%ld", edited.translations.size());
  draw_gllabel(win, "rotation", "%ld", edited.rotations.size());
  draw_gllabel(win, "scale", "%ld", edited.scales.size());
  draw_gllabel(win, "weights", "%ld", value.morphs.size());
  draw_gllabel(win, "targets", "%ld", edited.targets.size());
  if (updated) {
    edit = {sel.type, sel.index, edited, false};
  }
  return updated;
}

inline void draw_glscenetree(const opengl_window& win, const string& lbl,
    yocto_scene& scene, app_selection& sel, int height) {
  draw_glscenetree(win, scene, sel);
}

inline bool draw_glsceneinspector(const opengl_window& win, const string& lbl,
    yocto_scene& scene, app_selection& sel, app_edit& edit, int height) {
  if (sel.type == typeid(void)) return false;

  edit = app_edit{};
  if (sel.type == typeid(yocto_camera)) {
    return draw_glsceneinspector(
        win, scene.cameras[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_shape)) {
    return draw_glsceneinspector(
        win, scene.shapes[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_subdiv)) {
    return draw_glsceneinspector(
        win, scene.subdivs[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_texture)) {
    return draw_glsceneinspector(
        win, scene.textures[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_voltexture)) {
    return draw_glsceneinspector(
        win, scene.voltextures[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_material)) {
    return draw_glsceneinspector(
        win, scene.materials[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_environment)) {
    return draw_glsceneinspector(
        win, scene.environments[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_instance)) {
    return draw_glsceneinspector(
        win, scene.instances[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_scene_node)) {
    return draw_glsceneinspector(win, scene.nodes[sel.index], sel, edit, scene);
  }
  if (sel.type == typeid(yocto_animation)) {
    return draw_glsceneinspector(
        win, scene.animations[sel.index], sel, edit, scene);
  }

  return false;
}

#endif
