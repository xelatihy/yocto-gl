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

//
// =============================================================================
//
// WARNING: THIS IS REALLY UGLY CODE. THE DISNEY SCENE DATA HAS SPECIAL CASES
// ALMOST EVERYWHERE. THOSE SPECIAL CASES ARE HARDCODED IN THE CODE BELOW.
// REALLY, DO NOT USE THIS CODE.
//
// =============================================================================
//

#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
using namespace yocto;

#include "ext/CLI11.hpp"

#include "ext/json.hpp"
#include "ext/sajson.h"

#include <string_view>
#include <unordered_set>

using std::string_view;
using std::unordered_set;
using namespace std::literals::string_view_literals;

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

static inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)
    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get directory name (including '/').
static inline string get_dirname(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
static inline string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos + 1);
}

// Get filename without directory.
static inline string get_filename(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
static inline string get_noextension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
static inline string get_basename(const string& filename) {
  return get_noextension(get_filename(filename));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF DISNEY ISLAND SCENE
// -----------------------------------------------------------------------------
namespace yocto {

inline string replace(string_view str, string_view from, string_view to) {
  // https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
  auto replaced = ""s;
  while (!str.empty()) {
    auto pos = str.find(from);
    if (pos == string_view::npos) {
      replaced += str;
      break;
    } else {
      replaced += str.substr(0, pos);
      replaced += to;
      str.remove_prefix(pos + from.size());
    }
  }
  return replaced;
}

// Load a text file
static inline void load_text(const string& filename, string& str) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    throw std::runtime_error("cannot read file " + filename);
  }
  fclose(fs);
}

// Save a text file
inline void save_text(const string& filename, const string& str) {
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    throw std::runtime_error("cannot write file " + filename);
  }
  fclose(fs);
}

using json = nlohmann::json;

// Load a JSON object
inline void load_json(const string& filename, json& js) {
  auto text = ""s;
  load_text(filename, text);
  js = json::parse(text);
}

// Save a JSON object
inline void save_json(const string& filename, const json& js) {
  // we have to use streams here since the json library is faster with them
  save_text(filename, js.dump(4));
}

inline void to_json(json& js, const vec2f& val) {
  nlohmann::to_json(js, (const std::array<float, 2>&)val);
}
inline void from_json(const json& js, vec2f& val) {
  nlohmann::from_json(js, (std::array<float, 2>&)val);
}

inline void to_json(json& js, const vec3f& val) {
  nlohmann::to_json(js, (const std::array<float, 3>&)val);
}
inline void from_json(const json& js, vec3f& val) {
  nlohmann::from_json(js, (std::array<float, 3>&)val);
}

inline void to_json(json& js, const vec4f& val) {
  nlohmann::to_json(js, (const std::array<float, 4>&)val);
}
inline void from_json(const json& js, vec4f& val) {
  nlohmann::from_json(js, (std::array<float, 4>&)val);
}

inline void to_json(json& js, const vec2i& val) {
  nlohmann::to_json(js, (const std::array<int, 2>&)val);
}
inline void from_json(const json& js, vec2i& val) {
  nlohmann::from_json(js, (std::array<int, 2>&)val);
}

inline void to_json(json& js, const vec3i& val) {
  nlohmann::to_json(js, (const std::array<int, 3>&)val);
}
inline void from_json(const json& js, vec3i& val) {
  nlohmann::from_json(js, (std::array<int, 3>&)val);
}

inline void to_json(json& js, const vec4i& val) {
  nlohmann::to_json(js, (const std::array<int, 4>&)val);
}
inline void from_json(const json& js, vec4i& val) {
  nlohmann::from_json(js, (std::array<int, 4>&)val);
}

inline void to_json(json& js, const mat4f& val) {
  nlohmann::to_json(js, (const std::array<float, 16>&)val);
}
inline void from_json(const json& js, mat4f& val) {
  nlohmann::from_json(js, (std::array<float, 16>&)val);
}

struct disney_material {
  string name                = "";
  vec3f  color               = zero3f;
  string color_map           = ""s;
  float  refractive          = 0;
  int    color_ptex_faces    = 0;
  int    color_ptex_rowfaces = 0;
  int    color_ptex_colfaces = 0;
  int    color_ptex_tilesize = 0;
  string color_map_baked     = "";
};

void load_island_cameras(
    const string& filename, const string& dirname, yocto_scene& scene) {
  printf("%s\n", filename.c_str());
  auto js = json{};
  load_json(dirname + filename, js);
  auto camera  = yocto_camera{};
  camera.uri   = "cameras/" + get_basename(filename) + ".yaml";
  camera.lens  = js.at("focalLength").get<float>() * 0.001f;
  camera.focus = js.at("centerOfInterest").get<float>();
  // camera.aperture  = js.at("lensRadius").get<float>();
  camera.film.y = camera.film.x / js.at("ratio").get<float>();
  auto from     = js.at("eye").get<vec3f>();
  auto to       = js.at("look").get<vec3f>();
  auto up       = js.at("up").get<vec3f>();
  camera.frame  = lookat_frame(from, to, up);
  camera.focus  = length(from - to);
  scene.cameras.push_back(camera);
}

void load_island_lights(
    const string& filename, const string& dirname, yocto_scene& scene) {
  printf("%s\n", filename.c_str());
  auto js = json{};
  load_json(dirname + filename, js);
  for (auto& [name, ljs] : js.items()) {
    if (ljs.at("type") == "quad") {
      auto material     = yocto_material{};
      material.uri      = "materials/lights/" + name + ".yaml";
      material.emission = xyz(ljs.at("color").get<vec4f>()) *
                          pow(2.0f, ljs.at("exposure").get<float>());
      scene.materials.push_back(material);
      auto shape  = yocto_shape{};
      shape.uri   = "shapes/lights/" + name + ".ply";
      auto params = procshape_params{};
      params.type = make_shape_type::quad;
      params.scale =
          (ljs.at("width").get<float>() + ljs.at("height").get<float>());
      make_improc(shape.triangles, shape.quads, shape.positions, shape.normals,
          shape.texcoords, params);
      scene.shapes.push_back(shape);
      auto instance     = yocto_instance{};
      instance.uri      = "instances/lights/" + name + ".yaml";
      instance.frame    = frame3f(ljs.at("translationMatrix").get<mat4f>());
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = (int)scene.materials.size() - 1;
      scene.instances.push_back(instance);
    } else if (ljs.at("type") == "dome") {
      auto texture = yocto_texture{};
      texture.uri  = ljs.at("map");
      load_image(dirname + texture.uri, texture.hdr);
      scene.textures.push_back(texture);
      auto environment     = yocto_environment{};
      environment.uri      = "environments/lights/" + name + ".yaml";
      environment.emission = xyz(ljs.at("color").get<vec4f>()) *
                             pow(2.0f, ljs.at("exposure").get<float>());
      environment.emission_tex = (int)scene.textures.size() - 1;
      environment.frame = frame3f(ljs.at("translationMatrix").get<mat4f>());
      scene.environments.push_back(environment);
    } else {
      throw std::runtime_error("unknown light type");
    }
  }
}

void load_island_materials(const string& filename, const string& dirname,
    unordered_map<string, disney_material>& mmap) {
  auto tjs = json{};
  load_json(dirname + "textures/textures.json", tjs);
  auto js = json{};
  load_json(dirname + filename, js);
  for (auto& [mname, mjs] : js.items()) {
    auto material       = disney_material{};
    material.name       = mname;
    auto base           = mjs.at("baseColor").get<vector<float>>();
    material.refractive = mjs.value("refractive", 0.0f);
    material.color      = {
        powf(base[0], 2.2f), powf(base[1], 2.2f), pow(base[2], 2.2f)};
    material.color_map  = mjs.at("colorMap");
    mmap[material.name] = material;
    if (material.color_map != ""s) {
      for (auto jass : mjs.at("assignment")) {
        auto ass_material      = disney_material{};
        ass_material.name      = mname + "_" + jass.get<string>();
        ass_material.color_map = material.color_map + jass.get<string>() +
                                 ".ptx";
        ass_material.color_ptex_faces =
            tjs.at(ass_material.color_map).at("numFaces").get<int>();
        ass_material.color_ptex_rowfaces =
            tjs.at(ass_material.color_map).at("basedRowLength").get<int>();
        ass_material.color_ptex_colfaces =
            tjs.at(ass_material.color_map).at("basedColLength").get<int>();
        ass_material.color_ptex_tilesize =
            tjs.at(ass_material.color_map).at("basedTileSize").get<int>();
        ass_material.color_map_baked =
            tjs.at(ass_material.color_map).at("bakedFilename").get<string>();
        ass_material.color = xyz(
            tjs.at(ass_material.color_map).at("color").get<vec4f>());
        mmap[ass_material.name] = ass_material;
      }
    }
  }
}

struct load_island_shape_callbacks : obj_callbacks {
  vector<yocto_shape>&                    shapes;
  vector<yocto_material>&                 materials;
  unordered_map<string, vector<vec2i>>&   smap;
  unordered_map<string, disney_material>& mmap;
  unordered_map<string, int>&             tmap;
  yocto_scene&                            scene;
  const string&                           filename;
  const string&                           parent_name;

  load_island_shape_callbacks(vector<yocto_shape>& shapes,
      vector<yocto_material>&                      materials,
      unordered_map<string, vector<vec2i>>&        smap,
      unordered_map<string, disney_material>&      mmap,
      unordered_map<string, int>& tmap, yocto_scene& scene,
      const string& filename, const string& parent_name)
      : shapes{shapes}
      , materials{materials}
      , smap{smap}
      , mmap{mmap}
      , tmap{tmap}
      , scene{scene}
      , filename{filename}
      , parent_name{parent_name} {}

  // obj vertices
  std::deque<vec3f> opos  = std::deque<vec3f>();
  std::deque<vec3f> onorm = std::deque<vec3f>();

  // vertex maps
  unordered_map<int, int> pos_map  = unordered_map<int, int>();
  unordered_map<int, int> norm_map = unordered_map<int, int>();

  // last material and group name
  string                  gname      = ""s;
  string                  mname      = ""s;
  bool                    split_next = false;
  vector<disney_material> dmaterials = {};

  // Add  vertices to the current shape
  void add_fvverts(const vector<obj_vertex>& verts) {
    for (auto& vert : verts) {
      if (!vert.position) continue;
      auto pos_it = pos_map.find(vert.position);
      if (pos_it != pos_map.end()) continue;
      auto nverts = (int)shapes.back().positions.size();
      pos_map.insert(pos_it, {vert.position, nverts});
      shapes.back().positions.push_back(opos.at(vert.position - 1));
    }
    for (auto& vert : verts) {
      if (!vert.normal) continue;
      auto norm_it = norm_map.find(vert.normal);
      if (norm_it != norm_map.end()) continue;
      auto nverts = (int)shapes.back().normals.size();
      norm_map.insert(norm_it, {vert.normal, nverts});
      shapes.back().normals.push_back(onorm.at(vert.normal - 1));
    }
  }

  int add_texture(const string& mapname) {
    if (tmap.find(mapname) == tmap.end()) {
      scene.textures.push_back({});
      scene.textures.back().uri = mapname;
      tmap[mapname]             = scene.textures.size() - 1;
    }
    return tmap.at(mapname);
  }

  void split_shape() {
    if (!split_next) return;
    split_next     = false;
    auto dmaterial = mmap.at(mname);
    if (dmaterial.color_map != "") {
      try {
        dmaterial = mmap.at(mname + "_" + gname);
      } catch (std::out_of_range& e) {
        throw;
      }
    }
    if (!shapes.empty() && materials.back().uri == dmaterial.name) return;
    dmaterials.push_back(dmaterial);
    materials.push_back({});
    materials.back().uri = dmaterial.name;
    if (dmaterial.color_map != "") {
      materials.back().diffuse     = {1, 1, 1};
      materials.back().diffuse_tex = add_texture(dmaterial.color_map_baked);
    } else if (dmaterial.refractive == 0) {
      materials.back().diffuse = dmaterial.color;
    } else {
      materials.back().specular     = {1, 1, 1};
      materials.back().transmission = {1, 1, 1};
      materials.back().roughness    = 0;
      materials.back().thin         = false;
    }
    shapes.push_back(yocto_shape{});
    shapes.back().uri = "shapes/" + parent_name + "/" + get_basename(filename) +
                        "-" + std::to_string((int)shapes.size()) + ".ply";
    pos_map.clear();
    pos_map.reserve(1024 * 1024);
    norm_map.clear();
    norm_map.reserve(1024 * 1024);
  }

  void vert(const vec3f& v) override { opos.push_back(v); }
  void norm(const vec3f& v) override { onorm.push_back(v); }
  void texcoord(const vec2f& v) override {
    throw std::runtime_error("texture coord not supported");
  }
  void face(const vector<obj_vertex>& verts) override {
    split_shape();
    add_fvverts(verts);
    if (verts.size() == 4) {
      shapes.back().quadspos.push_back(
          {pos_map.at(verts[0].position), pos_map.at(verts[1].position),
              pos_map.at(verts[2].position), pos_map.at(verts[3].position)});
      shapes.back().quadsnorm.push_back(
          {norm_map.at(verts[0].normal), norm_map.at(verts[1].normal),
              norm_map.at(verts[2].normal), norm_map.at(verts[3].normal)});
    } else {
      for (auto i = 2; i < verts.size(); i++)
        shapes.back().quadspos.push_back(
            {pos_map.at(verts[0].position), pos_map.at(verts[i - 1].position),
                pos_map.at(verts[i].position), pos_map.at(verts[i].position)});
      for (auto i = 2; i < verts.size(); i++)
        shapes.back().quadsnorm.push_back(
            {norm_map.at(verts[0].normal), norm_map.at(verts[i - 1].normal),
                norm_map.at(verts[i].normal), norm_map.at(verts[i].normal)});
    }
    if (dmaterials.back().color_ptex_faces) {
      auto offset = (int)shapes.back().texcoords.size();
      auto face   = (int)shapes.back().quadstexcoord.size();
      face %= dmaterials.back().color_ptex_faces;
      auto face_i = face % dmaterials.back().color_ptex_rowfaces;
      auto face_j = face / dmaterials.back().color_ptex_rowfaces;
      if (verts.size() == 4) {
        auto du  = 1 / (float)dmaterials.back().color_ptex_rowfaces;
        auto dv  = 1 / (float)dmaterials.back().color_ptex_colfaces;
        auto dpu = 1 / (float)(dmaterials.back().color_ptex_rowfaces *
                               dmaterials.back().color_ptex_tilesize);
        auto dpv = 1 / (float)(dmaterials.back().color_ptex_colfaces *
                               dmaterials.back().color_ptex_tilesize);
        auto u0 = (face_i + 0) * du + dpu, v0 = (face_j + 0) * dv + dpv;
        auto u1 = (face_i + 1) * du - dpu, v1 = (face_j + 1) * dv - dpv;
        shapes.back().texcoords.push_back({u0, v0});
        shapes.back().texcoords.push_back({u1, v0});
        shapes.back().texcoords.push_back({u1, v1});
        shapes.back().texcoords.push_back({u0, v1});
        shapes.back().quadstexcoord.push_back(
            {offset + 0, offset + 1, offset + 2, offset + 3});
      } else {
        throw std::runtime_error("BAD PTEX TEXCOORDS");
      }
    }
  }
  void group(const string& name) override {
    gname      = name;
    split_next = true;
  }
  void usemtl(const string& name) override {
    mname      = name;
    split_next = true;
  }
};

void add_island_shape(yocto_scene& scene, const string& parent_name,
    const string& filename, const string& dirname,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  if (smap.find(filename) != smap.end()) return;
  printf("%s\n", filename.c_str());

  auto shapes      = vector<yocto_shape>{};
  auto materials   = vector<yocto_material>{};
  auto facevarying = false;

  try {
    // load obj
    auto cb = load_island_shape_callbacks{
        shapes, materials, smap, mmap, tmap, scene, filename, parent_name};
    load_obj(dirname + filename, cb, true);

    // check for PTEX errors
    for (auto id = 0; id < shapes.size(); id++) {
      if (cb.dmaterials[id].color_map != "") {
        auto ptex_faces  = cb.dmaterials[id].color_ptex_faces;
        auto shape_faces = max(
            (int)shapes[id].quads.size(), (int)shapes[id].quadspos.size());
        auto is_multiple = shape_faces % ptex_faces == 0;
        if (!is_multiple)
          printf("PTEX ERROR:  %d %d\n", ptex_faces, shape_faces);
      }
    }

    // conversion to non-facevarying
    if (!facevarying) {
      for (auto& shape : shapes) {
        auto split_quads     = vector<vec4i>{};
        auto split_positions = vector<vec3f>{};
        auto split_normals   = vector<vec3f>{};
        auto split_texcoords = vector<vec2f>{};
        split_facevarying(split_quads, split_positions, split_normals,
            split_texcoords, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords);
        shape.quads         = split_quads;
        shape.positions     = split_positions;
        shape.normals       = split_normals;
        shape.texcoords     = split_texcoords;
        shape.quadspos      = {};
        shape.quadsnorm     = {};
        shape.quadstexcoord = {};
        if (shape.texcoords.empty()) {
          auto all_triangles = true;
          for (auto& q : shape.quads) {
            if (q.z != q.w) {
              all_triangles = false;
              break;
            }
          }
          if (all_triangles) {
            quads_to_triangles(shape.triangles, shape.quads);
            shape.quads = {};
          }
        }
      }
    } else {
      for (auto& shape : shapes) {
        shape.uri = get_noextension(shape.uri) + ".obj";
      }
    }

// merging quads and triangles
#if 0
        if(!facevarying) {
            for (auto& shape : shapes) {
                merge_triangles_and_quads(shape.triangles, shape.quads, false);
                if (shape.triangles.empty() && shape.quads.empty())
                    throw std::runtime_error("empty shape");
            }
        }
#endif

  } catch (const std::exception& e) {
    throw std::runtime_error("cannot load mesh " + filename + "\n" + e.what());
  }

  for (auto shape_id = 0; shape_id < shapes.size(); shape_id++) {
    scene.shapes.push_back(shapes[shape_id]);
    scene.materials.push_back(materials[shape_id]);
    smap[filename].push_back(
        {(int)scene.shapes.size() - 1, (int)scene.materials.size() - 1});
  }
}

void add_island_instance(yocto_scene& scene, const string& parent_name,
    const mat4f& xform, const vector<vec2i>& shapes) {
  static auto name_counter = unordered_map<string, int>{};
  for (auto shape_material : shapes) {
    auto instance = yocto_instance{};
    instance.uri  = "instances/" + parent_name + "/" + parent_name + "-" +
                   std::to_string(name_counter[parent_name]++) + ".yaml";
    instance.frame    = frame3f(xform);
    instance.shape    = shape_material.x;
    instance.material = shape_material.y;
    scene.instances.push_back(instance);
  }
}

void add_island_variant_instance(vector<yocto_instance>& instances,
    const string& parent_name, const mat4f& xform,
    const vector<vec2i>& shapes) {
  for (auto shape_material : shapes) {
    auto instance     = yocto_instance{};
    instance.frame    = frame3f(xform);
    instance.shape    = shape_material.x;
    instance.material = shape_material.y;
    instances.push_back(instance);
  }
}

void load_island_archive(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  printf("%s\n", filename.c_str());
  auto buffer = ""s;
  load_text(dirname + filename, buffer);
  auto view = sajson::mutable_string_view(buffer.size(), buffer.data());
  auto doc  = sajson::parse(sajson::dynamic_allocation(), view);
  auto iijs = doc.get_root();
  for (auto j = 0; j < iijs.get_length(); j++) {
    auto shape_filename = iijs.get_object_key(j).as_string();
    add_island_shape(
        scene, parent_name, shape_filename, dirname, smap, mmap, tmap);
    auto xforms = iijs.get_object_value(j);
    for (auto i = 0; i < xforms.get_length(); i++) {
      auto xform_ = xforms.get_object_value(i);
      auto xform  = mat4f{};
      for (auto c = 0; c < 16; c++) {
        (&xform.x.x)[c] = xform_.get_array_element(c).get_number_value();
      }
      xform = parent_xform * xform;
      add_island_instance(scene, parent_name, xform, smap.at(shape_filename));
    }
  }
}

void load_island_variant_archive(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    vector<yocto_instance>&                 instances,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  // elements
  printf("%s\n", filename.c_str());
  auto buffer = ""s;
  load_text(dirname + filename, buffer);
  auto view = sajson::mutable_string_view(buffer.size(), buffer.data());
  auto doc  = sajson::parse(sajson::dynamic_allocation(), view);
  auto iijs = doc.get_root();
  for (auto j = 0; j < iijs.get_length(); j++) {
    auto shape_filename = iijs.get_object_key(j).as_string();
    add_island_shape(
        scene, parent_name, shape_filename, dirname, smap, mmap, tmap);
    auto xforms = iijs.get_object_value(j);
    for (auto i = 0; i < xforms.get_length(); i++) {
      auto xform_ = xforms.get_object_value(i);
      auto xform  = mat4f{};
      for (auto c = 0; c < 16; c++) {
        (&xform.x.x)[c] = xform_.get_array_element(c).get_number_value();
      }
      xform = parent_xform * xform;
      add_island_variant_instance(
          instances, parent_name, xform, smap.at(shape_filename));
    }
  }
}

void load_island_variants(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    unordered_map<string, vector<yocto_instance>>& instances,
    unordered_map<string, vector<vec2i>>&          smap,
    unordered_map<string, disney_material>&        mmap,
    unordered_map<string, int>&                    tmap) {
  printf("%s\n", filename.c_str());
  auto js_ = json{};
  load_json(dirname + filename, js_);

  for (auto& [vname, vjs] : js_.at("variants").items()) {
    vjs["transformMatrix"] = js_.at("transformMatrix");
  }

  // main instance
  for (auto& [vname, vjs] : js_.at("variants").items()) {
    instances[vname] = {};
    add_island_shape(
        scene, vname, vjs.at("geomObjFile"), dirname, smap, mmap, tmap);
    add_island_variant_instance(instances[vname], vname,
        vjs.at("transformMatrix"), smap.at(vjs.at("geomObjFile")));

    // instanced archives
    for (auto& [iiname, ijs] : vjs.at("instancedPrimitiveJsonFiles").items()) {
      auto filename = ijs.at("jsonFile").get<string>();
      if (ijs.at("type") == "archive") {
        load_island_variant_archive(filename, dirname, scene, vname,
            vjs.at("transformMatrix"), instances[vname], smap, mmap, tmap);
      } else {
        throw std::runtime_error("unknown instance type");
      }
    }
  }
}

void load_island_element(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  unordered_map<string, vector<yocto_instance>> variants;
  load_island_variants("json/isBayCedarA1/isBayCedarA1.json", dirname, scene,
      parent_name, identity4x4f, variants, smap, mmap, tmap);

  printf("%s\n", filename.c_str());
  auto buffer = ""s;
  load_text(dirname + filename, buffer);
  auto view = sajson::mutable_string_view(buffer.size(), buffer.data());
  auto doc  = sajson::parse(sajson::dynamic_allocation(), view);
  auto iijs = doc.get_root();
  for (auto j = 0; j < iijs.get_length(); j++) {
    auto  vname   = iijs.get_object_key(j).as_string();
    auto& variant = variants.at(vname);
    auto  xforms  = iijs.get_object_value(j);
    for (auto i = 0; i < xforms.get_length(); i++) {
      auto xform_ = xforms.get_object_value(i);
      auto xform  = mat4f{};
      for (auto c = 0; c < 16; c++) {
        (&xform.x.x)[c] = xform_.get_array_element(c).get_number_value();
      }
      xform = parent_xform * xform;
      for (auto& instance : variant) {
        add_island_instance(scene, parent_name, xform * (mat4f)instance.frame,
            {{instance.shape, instance.material}});
      }
    }
  }
}

void load_island_curve(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    float start_radius, float end_radius,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  printf("%s\n", filename.c_str());
  auto buffer = ""s;
  load_text(dirname + filename, buffer);
  auto view    = sajson::mutable_string_view(buffer.size(), buffer.data());
  auto doc     = sajson::parse(sajson::dynamic_allocation(), view);
  auto outname = "ply/" + get_dirname(filename).substr(5) +
                 get_filename(filename) + ".ply";
  if (smap.find(outname) == smap.end()) {
    auto curves   = doc.get_root();
    auto shape    = yocto_shape{};
    auto material = -1;
    for (auto j = 0; j < curves.get_length(); j++) {
      auto curve = curves.get_array_element(j);
      shape.uri  = outname;
      for (auto i = 0; i < curve.get_length(); i++) {
        auto point = curve.get_array_element(i);
        shape.positions.push_back({
            (float)point.get_array_element(0).get_number_value(),
            (float)point.get_array_element(1).get_number_value(),
            (float)point.get_array_element(2).get_number_value(),
        });
        shape.radius.push_back(
            lerp(start_radius, end_radius, (float)j / curve.get_length()));
        if (i != 0) {
          shape.lines.push_back({(int)shape.positions.size() - 2,
              (int)shape.positions.size() - 1});
        }
      }
    }
    scene.shapes.push_back(shape);
    smap[outname] = {{(int)scene.shapes.size() - 1, material}};
  }
  add_island_instance(scene, parent_name, parent_xform, smap.at(outname));
}

void load_island_curvetube(const string& filename, const string& dirname,
    yocto_scene& scene, const string& parent_name, const mat4f& parent_xform,
    float start_width, float end_width, const string& material_name,
    unordered_map<string, vector<vec2i>>&   smap,
    unordered_map<string, disney_material>& mmap,
    unordered_map<string, int>&             tmap) {
  printf("%s\n", filename.c_str());
  auto buffer = ""s;
  load_text(dirname + filename, buffer);
  auto view    = sajson::mutable_string_view(buffer.size(), buffer.data());
  auto doc     = sajson::parse(sajson::dynamic_allocation(), view);
  auto outname = "ply/" + get_dirname(filename).substr(5) +
                 get_filename(filename) + ".ply";
  if (smap.find(outname) == smap.end()) {
    auto curves      = doc.get_root();
    auto material    = yocto_material{};
    material.diffuse = mmap.at(material_name).color;
    scene.materials.push_back(material);
    auto shape      = yocto_shape{};
    shape.uri       = outname;
    auto ssmaterial = (int)scene.materials.size() - 1;
    for (auto j = 0; j < curves.get_length(); j++) {
      auto curve             = curves.get_array_element(j);
      auto bspline_positions = vector<vec3f>{};
      for (auto i = 0; i < curve.get_length(); i++) {
        auto point = curve.get_array_element(i);
        bspline_positions.push_back({
            (float)point.get_array_element(0).get_number_value(),
            (float)point.get_array_element(1).get_number_value(),
            (float)point.get_array_element(2).get_number_value(),
        });
        if (i == 0) {
          bspline_positions.push_back(bspline_positions.front());
        }
      }
      bspline_positions.push_back(bspline_positions.back());
      // bspline to cubic bezier from pbrt
      auto bezier_positions = vector<vec3f>{};
      for (auto i = 0; i < bspline_positions.size() - 3; i++) {
        // First compute equivalent Bezier control points.
        auto p01 = bspline_positions[i + 0];
        auto p12 = bspline_positions[i + 1];
        auto p23 = bspline_positions[i + 2];

        // We already have p12.
        auto p11 = lerp(p01, p12, 0.5f);
        auto p22 = lerp(p12, p23, 0.5f);

        // Now elevate to degree 3.
        if (i == 0) bezier_positions.push_back(p11);
        bezier_positions.push_back(lerp(p11, p12, 2 / 3.f));
        bezier_positions.push_back(lerp(p12, p22, 1 / 3.f));
        bezier_positions.push_back(p22);
      }
      auto bezier_radius = vector<float>{};
      for (auto i = 0; i < bezier_positions.size(); i++) {
        bezier_radius.push_back(lerp(start_width, end_width,
                                    (float)i / (bezier_positions.size() - 1)) /
                                2);
      }
      for (auto i = 0; i < (int)bezier_positions.size() - 1; i++) {
        auto p0 = bezier_positions[i], p1 = bezier_positions[i + 1];
        auto r0 = bezier_radius[i], r1 = bezier_radius[i + 1];
        auto h          = length(p1 - p0);
        auto f          = frame_fromz(p0, p1 - p0);
        auto qpositions = vector<vec3f>{{r0, 0, 0}, {0, r0, 0}, {-r0, 0, 0},
            {0, -r0, 0}, {r1, 0, h}, {0, r1, h}, {-r1, 0, h}, {0, -r1, h}};
        auto qnormals   = vector<vec3f>{{1, 0, 0}, {0, 1, 0}, {-1, 0, 0},
            {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}};
        for (auto& p : qpositions) p = transform_point(f, p);
        for (auto& n : qnormals) n = transform_direction(f, n);
        auto qquads = vector<vec4i>{
            {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
        merge_quads(shape.quads, shape.positions, shape.normals,
            shape.texcoords, qquads, qpositions, qnormals, {});
      }
    }
    scene.shapes.push_back(shape);
    smap[outname] = {{(int)scene.shapes.size() - 1, ssmaterial}};
  }
  add_island_instance(scene, parent_name, parent_xform, smap.at(outname));
}

void load_island_elements(const string& filename, const string& dirname,
    yocto_scene& scene, unordered_map<string, vector<vec2i>>& smap,
    unordered_map<string, int>& tmap) {
  // instancing model
  // - main shape: "geomObjFile" and "name" properties
  // - main instance: "transform"
  // - main instances: "instancedPrimitiveJsonFiles"
  // - instances can be "archive" -> list of instances
  //                    "curve"   -> list of curve data (one shape only)
  //                    "element" -> ??? (happens only once)
  // - instanced copies:
  //                    copy element json data
  //                    set the main "geomObjFile" if not present
  //                    check what to do with "instancedPrimitiveJsonFiles"
  // - shapes: to avoid duplication, create a map of shapes from obj paths
  // - variants: what are they?
  // - materials: material names are not absolute; prepend element name

  printf("%s\n", filename.c_str());
  auto js = json{};
  load_json(dirname + filename, js);

  // add empty elements for simplicity
  if (!js.count("instancedPrimitiveJsonFiles")) {
    js["instancedPrimitiveJsonFiles"] = json::object();
  }
  if (!js.count("instancedCopies")) {
    js["instancedCopies"] = json::object();
  }

  // materials
  auto mmap = unordered_map<string, disney_material>{};
  load_island_materials(js.at("matFile"), dirname, mmap);

  // main instance
  auto name = js.at("name").get<string>();
  add_island_shape(
      scene, name, js.at("geomObjFile"), dirname, smap, mmap, tmap);
  add_island_instance(
      scene, name, js.at("transformMatrix"), smap.at(js.at("geomObjFile")));

  // instanced archives
  for (auto& [iiname, ijs] : js.at("instancedPrimitiveJsonFiles").items()) {
    auto filename = ijs.at("jsonFile").get<string>();
    if (ijs.at("type") == "archive") {
      load_island_archive(filename, dirname, scene, name,
          js.at("transformMatrix"), smap, mmap, tmap);
    } else if (ijs.at("type") == "curve") {
      load_island_curvetube(filename, dirname, scene, name,
          js.at("transformMatrix"), ijs.at("widthRoot"), ijs.at("widthTip"),
          ijs.at("material"), smap, mmap, tmap);
    } else if (ijs.at("type") == "element") {
      load_island_element(filename, dirname, scene, name,
          js.at("transformMatrix"), smap, mmap, tmap);
    } else if (ijs.at("type") == "skip") {
      printf("skipping %s\n", filename.c_str());
    } else {
      throw std::runtime_error("unknown instance type");
    }
  }

  // instanced copies
  for (auto& [iname, cjs] : js.at("instancedCopies").items()) {
    if (cjs.count("geomObjFile")) {
      add_island_shape(
          scene, name, cjs.at("geomObjFile"), dirname, smap, mmap, tmap);
      add_island_instance(scene, name, cjs.at("transformMatrix"),
          smap.at(cjs.at("geomObjFile")));
    } else {
      add_island_instance(scene, name, cjs.at("transformMatrix"),
          smap.at(js.at("geomObjFile")));
    }
    if (cjs.count("instancedPrimitiveJsonFiles") ||
        js.count("instancedPrimitiveJsonFiles")) {
      for (auto& [iiname, ijs] :
          cjs.count("instancedPrimitiveJsonFiles")
              ? cjs.at("instancedPrimitiveJsonFiles").items()
              : js.at("instancedPrimitiveJsonFiles").items()) {
        auto filename = ijs.at("jsonFile").get<string>();
        if (ijs.at("type") == "archive") {
          load_island_archive(filename, dirname, scene, name,
              cjs.at("transformMatrix"), smap, mmap, tmap);
        } else if (ijs.at("type") == "curve") {
          load_island_curvetube(filename, dirname, scene, name,
              cjs.at("transformMatrix"), ijs.at("widthRoot"),
              ijs.at("widthTip"), ijs.at("material"), smap, mmap, tmap);
        } else if (ijs.at("type") == "element") {
        } else if (ijs.at("type") == "skip") {
          printf("skipping %s\n", filename.c_str());
        } else {
          throw std::runtime_error("unknown instance type");
        }
      }
    }
  }

  // rename materials and shapes
}

void load_textures(
    yocto_scene& scene, const string& dirname, const load_params& params);

void load_island_scene(
    const string& filename, yocto_scene& scene, const load_params& params) {
  try {
    auto js = json{};
    load_json(filename, js);
    auto dirname = get_dirname(filename);

    for (auto filename : js.at("cameras").get<vector<string>>()) {
      load_island_cameras(filename, dirname, scene);
    }
    auto smap = std::unordered_map<string, vector<vec2i>>{};
    auto tmap = std::unordered_map<string, int>{};
    for (auto filename : js.at("elements").get<vector<string>>()) {
      load_island_elements(filename, dirname, scene, smap, tmap);
    }
    for (auto filename : js.at("lights").get<vector<string>>()) {
      load_island_lights(filename, dirname, scene);
    }

    // load meshes and textures
    load_textures(scene, dirname, params);
  } catch (std::exception& e) {
    throw std::runtime_error("error loading scene "s + e.what());
  }

  // fix texture names
  for (auto& texture : scene.textures) {
    texture.uri = replace(texture.uri, "ptex2png/", "textures/");
    texture.uri = replace(texture.uri, ".exr", ".hdr");
  }

  // fix names
  auto parent_shape_map = unordered_map<string, vec2i>{};
  for (auto id = 0; id < scene.shapes.size(); id++) {
    auto parent_name = get_dirname(scene.shapes[id].uri).substr(7);
    parent_name      = parent_name.substr(0, parent_name.size() - 1);
    parent_shape_map[parent_name].y += 1;
  }
  for (auto id = 0; id < scene.shapes.size(); id++) {
    auto parent_name = get_dirname(scene.shapes[id].uri).substr(7);
    parent_name      = parent_name.substr(0, parent_name.size() - 1);
    if (parent_shape_map[parent_name].y == 1) {
      scene.shapes[id].uri    = "shapes/" + parent_name + ".ply";
      scene.materials[id].uri = "materials/" + parent_name + ".yaml";
    } else {
      scene.shapes[id].uri = "shapes/" + parent_name +
                             std::to_string(parent_shape_map[parent_name].x) +
                             ".ply";
      scene.materials[id].uri =
          "materials/" + parent_name +
          std::to_string(parent_shape_map[parent_name].x) + ".ply";
      parent_shape_map[parent_name].x += 1;
    }
  }
  auto parent_texture_map = unordered_map<string, vec2i>{};
  for (auto id = 0; id < scene.textures.size(); id++) {
    auto parent_name = get_dirname(scene.textures[id].uri).substr(9);
    parent_name      = parent_name.substr(0, parent_name.size() - 1);
    parent_name      = replace(parent_name, "/Color", "");
    parent_texture_map[parent_name].y += 1;
  }
  for (auto id = 0; id < scene.textures.size(); id++) {
    auto parent_name = get_dirname(scene.textures[id].uri).substr(9);
    parent_name      = parent_name.substr(0, parent_name.size() - 1);
    if (parent_name == "lights") {
      scene.textures[id].uri = "textures/" +
                               get_filename(scene.textures[id].uri);
      continue;
    }
    parent_name = replace(parent_name, "/Color", "");
    if (parent_texture_map[parent_name].y == 1) {
      scene.textures[id].uri = "textures/" + parent_name + ".png";
    } else {
      scene.textures[id].uri =
          "textures/" + parent_name +
          std::to_string(parent_texture_map[parent_name].x) + ".png";
      parent_texture_map[parent_name].x += 1;
    }
  }
  // fix instances
  rename_instances(scene);

  // fix scene
  if (scene.uri == "") scene.uri = get_filename(filename);
  add_cameras(scene);
  add_materials(scene);
  normalize_uris(scene);
  trim_memory(scene);
  update_transforms(scene);

  // print stats
  printf("%s\n", format_stats(scene).c_str());
}

}  // namespace yocto

bool mkdir(const string& dir) {
  if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
    return true;
#ifndef _MSC_VER
  system(("mkdir -p " + dir).c_str());
  return true;
#else
  system(("mkdir " + dir).c_str());
  return true;
#endif
}

int main(int argc, char** argv) {
  // command line parameters
  auto notextures     = false;
  auto mesh_filenames = true;
  auto mesh_directory = "shapes/"s;
  auto uniform_txt    = false;
  auto validate       = false;
  auto info           = false;
  auto output         = "out.json"s;
  auto filename       = "scene.json"s;

  // parse command line
  auto parser = CLI::App{"Process scene"};
  parser.add_flag("--notextures", notextures, "Disable textures.");
  parser.add_flag("--mesh-filenames,!--no-mesh-filenames", mesh_filenames,
      "Add mesh filenames.");
  parser.add_option(
      "--mesh-directory", mesh_directory, "Mesh directory when adding names.");
  parser.add_flag("--uniform-textures,!--no-uniform-textures", uniform_txt,
      "uniform texture formats");
  parser.add_flag("--info,-i", info, "print scene info");
  parser.add_flag("--validate,!--no-validate", validate, "Validate scene");
  parser.add_option("--output,-o", output, "output scene")->required(true);
  parser.add_option("scene", filename, "input scene")->required(true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }
  setbuf(stdout, nullptr);

  // fix params
  auto load_prms       = load_params();
  auto save_prms       = save_params();
  load_prms.notextures = notextures;
  save_prms.notextures = notextures;

  // load scene
  auto scene = yocto_scene{};
  printf("loading scene");
  auto load_timer = timer();
  try {
    load_island_scene(filename, scene, load_prms);
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }
  printf(" in %s\n", load_timer.elapsedf().c_str());

  // validate scene
  if (validate) {
    printf("validating scene");
    auto validate_timer = timer();
    print_validation(scene);
    printf(" in %s\n", validate_timer.elapsedf().c_str());
  }

  // print info
  if (info) printf("%s\n", format_stats(scene).c_str());

// add missing mesh names if necessary
#if 0
    if (!mesh_directory.empty() && mesh_directory.back() != '/')
        mesh_directory += '/';
    if (mesh_filenames && get_extension(output) == "json") {
        for (auto& shape : scene.shapes) {
            shape.name = "";
            if (shape.positions.size() <= 16) continue;
            if (shape.facevarying) {
                shape.filename = mesh_directory + shape.name + ".obj";
            } else {
                shape.filename = mesh_directory + shape.name + ".ply";
            }
        }
    }
#endif

  // make a directory if needed
  auto dirname  = get_dirname(output);
  auto dirnames = unordered_set<string>{dirname};
  for (auto& shape : scene.shapes)
    dirnames.insert(dirname + get_dirname(shape.uri));
  for (auto& texture : scene.textures)
    dirnames.insert(dirname + get_dirname(texture.uri));
  if (get_extension(output) == "yaml") dirnames.insert(dirname + "instances/");
  for (auto& dir : dirnames) {
    if (!mkdir(get_dirname(dir))) {
      printf("cannot create directory %s\n", get_dirname(output).c_str());
      exit(1);
    }
  }

  // save scene
  printf("saving scene");
  auto save_timer = timer();
  try {
    save_prms.notextures = false;
    save_prms.noparallel = false;
    // save_prms.ply_instances = true;
    save_scene(output, scene, save_prms);
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }
  printf(" in %s\n", save_timer.elapsedf().c_str());

  // done
  return 0;
}
