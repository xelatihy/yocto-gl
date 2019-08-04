//
// Implementation for Yocto/Obj.
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_obj.h"

#include <algorithm>
#include <string_view>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;
using namespace std::literals::string_view_literals;

// Check if a file can be opened for reading.
static inline bool exists_file(const string& filename) {
  auto f = fopen(filename.c_str(), "r");
  if (!f) return false;
  fclose(f);
  return true;
}

// A file holder that closes a file when destructed. Useful for RIIA
struct file_holder {
  FILE*  fs       = nullptr;
  string filename = "";

  file_holder(const file_holder&) = delete;
  file_holder& operator=(const file_holder&) = delete;
  ~file_holder() {
    if (fs) fclose(fs);
  }
};

// Opens a file returing a handle with RIIA
static inline file_holder open_input_file(
    const string& filename, bool binary = false) {
  auto fs = fopen(filename.c_str(), !binary ? "rt" : "rb");
  if (!fs) throw std::runtime_error("could not open file " + filename);
  return {fs, filename};
}

// Read a line
static inline bool read_line(FILE* fs, char* buffer, size_t size) {
  return fgets(buffer, size, fs) != nullptr;
}

static inline bool is_obj_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_obj_newline(char c) { return c == '\r' || c == '\n'; }

static inline void skip_obj_whitespace(string_view& str) {
  while (!str.empty() && is_obj_space(str.front())) str.remove_prefix(1);
}
static inline void remove_obj_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_obj_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
static inline void parse_obj_value(string_view& str, string_view& value) {
  skip_obj_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_obj_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
    if (str.front() != '"') throw std::runtime_error("cannot parse value");
    str.remove_prefix(1);
    if (str.empty()) throw std::runtime_error("cannot parse value");
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
    if (cpy.empty()) throw std::runtime_error("cannot parse value");
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    str.remove_prefix(1);
  }
}
static inline void parse_obj_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_obj_value(str, valuev);
  value = string{valuev};
}
static inline void parse_obj_value(string_view& str, int& value) {
  char* end = nullptr;
  value     = (int)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_obj_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_obj_value(str, valuei);
  value = (bool)valuei;
}
static inline void parse_obj_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
static inline void parse_obj_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_obj_value(str, values[i]);
}

static inline void parse_obj_value(string_view& str, vec2f& value) {
  parse_obj_value(str, &value.x, 2);
}
static inline void parse_obj_value(string_view& str, vec3f& value) {
  parse_obj_value(str, &value.x, 3);
}
static inline void parse_obj_value(string_view& str, frame3f& value) {
  parse_obj_value(str, &value.x.x, 12);
}

template <typename T>
static inline void parse_obj_value_or_empty(string_view& str, T& value) {
  skip_obj_whitespace(str);
  if (str.empty()) {
    value = T{};
  } else {
    parse_obj_value(str, value);
  }
}

static inline void parse_obj_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  parse_obj_value(str, value.position);
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      parse_obj_value(str, value.normal);
    } else {
      parse_obj_value(str, value.texcoord);
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        parse_obj_value(str, value.normal);
      }
    }
  }
}

// Input for OBJ textures
static inline void parse_obj_value(string_view& str, mtl_texture_info& info) {
  // initialize
  info = mtl_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_obj_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    parse_obj_value(str, token);
    tokens.push_back(token);
    skip_obj_whitespace(str);
  }
  if (tokens.empty()) throw std::runtime_error("cannot parse value");

  // texture name
  info.path = fs::path(tokens.back()).generic_string();

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }
}

// Load obj materials
void load_mtl(const string& filename, obj_callbacks& cb, bool fliptr) {
  // open file
  auto fs_ = open_input_file(filename);
  auto fs  = fs_.fs;

  // currently parsed material
  auto material = mtl_material{};
  auto first    = true;

  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      if (!first) cb.material(material);
      first    = false;
      material = mtl_material{};
      parse_obj_value(line, material.name);
    } else if (cmd == "illum") {
      parse_obj_value(line, material.illum);
    } else if (cmd == "Ke") {
      parse_obj_value(line, material.ke);
    } else if (cmd == "Kd") {
      parse_obj_value(line, material.kd);
    } else if (cmd == "Ks") {
      parse_obj_value(line, material.ks);
    } else if (cmd == "Kt") {
      parse_obj_value(line, material.kt);
    } else if (cmd == "Tf") {
      material.kt = {-1, -1, -1};
      parse_obj_value(line, material.kt);
      if (material.kt.y < 0)
        material.kt = {material.kt.x, material.kt.x, material.kt.x};
      if (fliptr) material.kt = vec3f{1, 1, 1} - material.kt;
    } else if (cmd == "Tr") {
      parse_obj_value(line, material.op);
      if (fliptr) material.op = 1 - material.op;
    } else if (cmd == "Ns") {
      parse_obj_value(line, material.ns);
      material.pr = pow(2 / (material.ns + 2), 1 / 4.0f);
      if (material.pr < 0.01f) material.pr = 0;
      if (material.pr > 0.99f) material.pr = 1;
    } else if (cmd == "d") {
      parse_obj_value(line, material.op);
    } else if (cmd == "map_Ke") {
      parse_obj_value(line, material.ke_map);
    } else if (cmd == "map_Kd") {
      parse_obj_value(line, material.kd_map);
    } else if (cmd == "map_Ks") {
      parse_obj_value(line, material.ks_map);
    } else if (cmd == "map_Tr") {
      parse_obj_value(line, material.kt_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_obj_value(line, material.op_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_obj_value(line, material.bump_map);
    } else if (cmd == "map_occ" || cmd == "occ") {
      parse_obj_value(line, material.occ_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_obj_value(line, material.disp_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_obj_value(line, material.norm_map);
    } else if (cmd == "Pm") {
      parse_obj_value(line, material.pm);
    } else if (cmd == "Pr") {
      parse_obj_value(line, material.pr);
    } else if (cmd == "Ps") {
      parse_obj_value(line, material.ps);
    } else if (cmd == "Pc") {
      parse_obj_value(line, material.pc);
    } else if (cmd == "Pcr") {
      parse_obj_value(line, material.pcr);
    } else if (cmd == "map_Pm") {
      parse_obj_value(line, material.pm_map);
    } else if (cmd == "map_Pr") {
      parse_obj_value(line, material.pr_map);
    } else if (cmd == "map_Ps") {
      parse_obj_value(line, material.ps_map);
    } else if (cmd == "map_Pc") {
      parse_obj_value(line, material.pc_map);
    } else if (cmd == "map_Pcr") {
      parse_obj_value(line, material.pcr_map);
    } else if (cmd == "Vt") {
      parse_obj_value(line, material.vt);
    } else if (cmd == "Vp") {
      parse_obj_value(line, material.vp);
    } else if (cmd == "Ve") {
      parse_obj_value(line, material.ve);
    } else if (cmd == "Vs") {
      parse_obj_value(line, material.vs);
    } else if (cmd == "Vg") {
      parse_obj_value(line, material.vg);
    } else if (cmd == "Vr") {
      parse_obj_value(line, material.vr);
    } else if (cmd == "map_Vs") {
      parse_obj_value(line, material.vs_map);
    }
  }

  // issue current material
  if (!first) cb.material(material);
}

// Load obj extensions
void load_objx(const string& filename, obj_callbacks& cb) {
  // open file
  auto fs_ = open_input_file(filename);
  auto fs  = fs_.fs;

  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "c") {
      auto camera = objx_camera();
      parse_obj_value(line, camera.name);
      parse_obj_value(line, camera.ortho);
      parse_obj_value(line, camera.width);
      parse_obj_value(line, camera.height);
      parse_obj_value(line, camera.lens);
      parse_obj_value(line, camera.focus);
      parse_obj_value(line, camera.aperture);
      parse_obj_value(line, camera.frame);
      cb.camera(camera);
    } else if (cmd == "e") {
      auto environment = objx_environment();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.ke);
      parse_obj_value(line, environment.ke_txt.path);
      parse_obj_value(line, environment.frame);
      if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
      cb.environmnet(environment);
    } else if (cmd == "i") {
      auto instance = objx_instance();
      parse_obj_value(line, instance.name);
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.material);
      parse_obj_value(line, instance.frame);
      cb.instance(instance);
    } else if (cmd == "po") {
      auto procedural = objx_procedural();
      parse_obj_value(line, procedural.name);
      parse_obj_value(line, procedural.type);
      parse_obj_value(line, procedural.material);
      parse_obj_value(line, procedural.size);
      parse_obj_value(line, procedural.level);
      parse_obj_value(line, procedural.frame);
      cb.procedural(procedural);
    } else {
      // unused
    }
  }
}

// Load obj scene
void load_obj(const string& filename, obj_callbacks& cb, bool nomaterials,
    bool flipv, bool fliptr) {
  // open file
  auto fs_ = open_input_file(filename);
  auto fs  = fs_.fs;

  // track vertex size
  auto vert_size = obj_vertex();
  auto verts     = vector<obj_vertex>();  // buffer to avoid reallocation

  // material libraries read already
  auto mlibs = vector<string>{};

  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      auto vert = zero3f;
      parse_obj_value(line, vert);
      cb.vert(vert);
      vert_size.position += 1;
    } else if (cmd == "vn") {
      auto vert = zero3f;
      parse_obj_value(line, vert);
      cb.norm(vert);
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      auto vert = zero2f;
      parse_obj_value(line, vert);
      if (flipv) vert.y = 1 - vert.y;
      cb.texcoord(vert);
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      verts.clear();
      skip_obj_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_obj_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        verts.push_back(vert);
        skip_obj_whitespace(line);
      }
      if (cmd == "f") cb.face(verts);
      if (cmd == "l") cb.line(verts);
      if (cmd == "p") cb.point(verts);
    } else if (cmd == "o") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.object(name);
    } else if (cmd == "usemtl") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.usemtl(name);
    } else if (cmd == "g") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.group(name);
    } else if (cmd == "s") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.smoothing(name);
    } else if (cmd == "mtllib") {
      if (nomaterials) continue;
      auto mtlname = ""s;
      parse_obj_value(line, mtlname);
      cb.mtllib(mtlname);
      if (std::find(mlibs.begin(), mlibs.end(), mtlname) != mlibs.end())
        continue;
      mlibs.push_back(mtlname);
      auto mtlpath = fs::path(filename).parent_path() / mtlname;
      load_mtl(mtlpath, cb, fliptr);
    } else {
      // unused
    }
  }

  // parse extensions if presents
  if (!nomaterials) {
    auto extname    = fs::path(filename).replace_extension(".objx");
    auto ext_exists = exists_file(extname);
    if (ext_exists) {
      load_objx(extname, cb);
    }
  }
}

// Read obj
bool read_obj_element(FILE* fs, obj_element& element, vec3f& value,
    string& name, vector<obj_vertex>& vertices, obj_vertex& vert_size) {
  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      element = obj_element::vertex;
      parse_obj_value(line, value);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      element = obj_element::normal;
      parse_obj_value(line, value);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      element = obj_element::texcoord;
      parse_obj_value(line, (vec2f&)value);
      vert_size.texcoord += 1;
      return true;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      skip_obj_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_obj_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        vertices.push_back(vert);
        skip_obj_whitespace(line);
      }
      if (cmd == "f") element = obj_element::face;
      if (cmd == "l") element = obj_element::line;
      if (cmd == "p") element = obj_element::point;
      return true;
    } else if (cmd == "o") {
      element = obj_element::object;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "usemtl") {
      element = obj_element::usemtl;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "g") {
      element = obj_element::group;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "s") {
      element = obj_element::smoothing;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "mtllib") {
      element = obj_element::mtllib;
      parse_obj_value(line, name);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
bool read_mtl_element(
    FILE* fs, mtl_element& element, mtl_material& material, bool fliptr) {
  // currently parsed material
  material   = mtl_material{};
  auto found = false;

  // read the file line by line
  auto fpos = ftell(fs);
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      if (found) {
        fseek(fs, fpos, SEEK_SET);
        element = mtl_element::material;
        return true;
      }
      found = true;
      parse_obj_value(line, material.name);
    } else if (cmd == "illum") {
      parse_obj_value(line, material.illum);
    } else if (cmd == "Ke") {
      parse_obj_value(line, material.ke);
    } else if (cmd == "Kd") {
      parse_obj_value(line, material.kd);
    } else if (cmd == "Ks") {
      parse_obj_value(line, material.ks);
    } else if (cmd == "Kt") {
      parse_obj_value(line, material.kt);
    } else if (cmd == "Tf") {
      material.kt = {-1, -1, -1};
      parse_obj_value(line, material.kt);
      if (material.kt.y < 0)
        material.kt = {material.kt.x, material.kt.x, material.kt.x};
      if (fliptr) material.kt = vec3f{1, 1, 1} - material.kt;
    } else if (cmd == "Tr") {
      parse_obj_value(line, material.op);
      if (fliptr) material.op = 1 - material.op;
    } else if (cmd == "Ns") {
      parse_obj_value(line, material.ns);
      material.pr = pow(2 / (material.ns + 2), 1 / 4.0f);
      if (material.pr < 0.01f) material.pr = 0;
      if (material.pr > 0.99f) material.pr = 1;
    } else if (cmd == "d") {
      parse_obj_value(line, material.op);
    } else if (cmd == "map_Ke") {
      parse_obj_value(line, material.ke_map);
    } else if (cmd == "map_Kd") {
      parse_obj_value(line, material.kd_map);
    } else if (cmd == "map_Ks") {
      parse_obj_value(line, material.ks_map);
    } else if (cmd == "map_Tr") {
      parse_obj_value(line, material.kt_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_obj_value(line, material.op_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_obj_value(line, material.bump_map);
    } else if (cmd == "map_occ" || cmd == "occ") {
      parse_obj_value(line, material.occ_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_obj_value(line, material.disp_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_obj_value(line, material.norm_map);
    } else if (cmd == "Pm") {
      parse_obj_value(line, material.pm);
    } else if (cmd == "Pr") {
      parse_obj_value(line, material.pr);
    } else if (cmd == "Ps") {
      parse_obj_value(line, material.ps);
    } else if (cmd == "Pc") {
      parse_obj_value(line, material.pc);
    } else if (cmd == "Pcr") {
      parse_obj_value(line, material.pcr);
    } else if (cmd == "map_Pm") {
      parse_obj_value(line, material.pm_map);
    } else if (cmd == "map_Pr") {
      parse_obj_value(line, material.pr_map);
    } else if (cmd == "map_Ps") {
      parse_obj_value(line, material.ps_map);
    } else if (cmd == "map_Pc") {
      parse_obj_value(line, material.pc_map);
    } else if (cmd == "map_Pcr") {
      parse_obj_value(line, material.pcr_map);
    } else if (cmd == "Vt") {
      parse_obj_value(line, material.vt);
    } else if (cmd == "Vp") {
      parse_obj_value(line, material.vp);
    } else if (cmd == "Ve") {
      parse_obj_value(line, material.ve);
    } else if (cmd == "Vs") {
      parse_obj_value(line, material.vs);
    } else if (cmd == "Vg") {
      parse_obj_value(line, material.vg);
    } else if (cmd == "Vr") {
      parse_obj_value(line, material.vr);
    } else if (cmd == "map_Vs") {
      parse_obj_value(line, material.vs_map);
    }

    // update pos
    fpos = ftell(fs);
  }

  // return found value
  if (found) {
    element = mtl_element::material;
    return true;
  } else {
    return false;
  }
}

// Read objx
bool read_objx_element(FILE* fs, objx_element& element, objx_camera& camera,
    objx_environment& environment, objx_instance& instance,
    objx_procedural& procedural) {
  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get element
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "c") {
      element = objx_element::camera;
      camera  = objx_camera();
      parse_obj_value(line, camera.name);
      parse_obj_value(line, camera.ortho);
      parse_obj_value(line, camera.width);
      parse_obj_value(line, camera.height);
      parse_obj_value(line, camera.lens);
      parse_obj_value(line, camera.focus);
      parse_obj_value(line, camera.aperture);
      parse_obj_value(line, camera.frame);
      return true;
    } else if (cmd == "e") {
      element     = objx_element::environment;
      environment = objx_environment();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.ke);
      parse_obj_value(line, environment.ke_txt.path);
      parse_obj_value(line, environment.frame);
      if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
      return true;
    } else if (cmd == "i") {
      element  = objx_element::instance;
      instance = objx_instance();
      parse_obj_value(line, instance.name);
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.material);
      parse_obj_value(line, instance.frame);
      return true;
    } else if (cmd == "po") {
      element    = objx_element::procedural;
      procedural = objx_procedural();
      parse_obj_value(line, procedural.name);
      parse_obj_value(line, procedural.type);
      parse_obj_value(line, procedural.material);
      parse_obj_value(line, procedural.size);
      parse_obj_value(line, procedural.level);
      parse_obj_value(line, procedural.frame);
      return true;
    } else {
      // unused
    }
  }

  return false;
}

// Write text to file
static inline void write_obj_value(FILE* fs, int value) {
  if (fprintf(fs, "%d", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, float value) {
  if (fprintf(fs, "%g", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_text(FILE* fs, const char* value) {
  if (fprintf(fs, "%s", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_text(FILE* fs, const string& value) {
  if (fprintf(fs, "%s", value.c_str()) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, const char* value) {
  if (fprintf(fs, "%s", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, const string& value) {
  if (fprintf(fs, "%s", value.c_str()) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, const vec2f& value) {
  if (fprintf(fs, "%g %g", value.x, value.y) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, const vec3f& value) {
  if (fprintf(fs, "%g %g %g", value.x, value.y, value.z) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_obj_value(FILE* fs, const frame3f& value) {
  for (auto i = 0; i < 12; i++)
    if (fprintf(fs, i ? " %g" : "%g", (&value.x.x)[i]) < 0)
      throw std::runtime_error("cannot print value");
}
static void write_obj_value(FILE* fs, const obj_vertex& value) {
  if (fprintf(fs, "%d", value.position) < 0)
    throw std::runtime_error("cannot write value");
  if (value.texcoord) {
    if (fprintf(fs, "/%d", value.texcoord) < 0)
      throw std::runtime_error("cannot write value");
    if (value.normal) {
      if (fprintf(fs, "/%d", value.normal) < 0)
        throw std::runtime_error("cannot write value");
    }
  } else if (value.normal) {
    if (fprintf(fs, "//%d", value.normal) < 0)
      throw std::runtime_error("cannot write value");
  }
}

template <typename T, typename... Ts>
static inline void write_obj_line_(
    FILE* fs, const T& value, const Ts... values) {
  write_obj_value(fs, value);
  if constexpr (sizeof...(values) == 0) {
    write_obj_text(fs, "\n");
  } else {
    write_obj_text(fs, " ");
    write_obj_line_(fs, values...);
  }
}
template <typename T, typename Ts>
static inline void write_obj_line_(
    FILE* fs, const T& value, const vector<Ts>& values) {
  write_obj_value(fs, value);
  for (auto& value : values) {
    write_obj_text(fs, " ");
    write_obj_value(fs, value);
  }
  write_obj_text(fs, "\n");
}

static inline vector<string> split_string(
    const string& str, const string& delim) {
  auto tokens = vector<string>{};
  auto last = (size_t)0, next = (size_t)0;
  while ((next = str.find(delim, last)) != string::npos) {
    tokens.push_back(str.substr(last, next - last));
    last = next + delim.size();
  }
  if (last < str.size()) tokens.push_back(str.substr(last));
  return tokens;
}

// Write obj elements
void write_obj_comment(FILE* fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    write_obj_text(fs, "# " + line + "\n");
  }
  write_obj_text(fs, "\n");
}

void write_obj_element(FILE* fs, obj_element element, const vec3f& value,
    const string& name, const vector<obj_vertex>& vertices) {
  switch (element) {
    case obj_element::vertex: write_obj_line_(fs, "v", value); break;
    case obj_element::normal: write_obj_line_(fs, "vn", value); break;
    case obj_element::texcoord:
      write_obj_line_(fs, "vt", vec2f{value.x, value.y});
      break;
    case obj_element::face: write_obj_line_(fs, "f", vertices); break;
    case obj_element::line: write_obj_line_(fs, "l", vertices); break;
    case obj_element::point: write_obj_line_(fs, "p", vertices); break;
    case obj_element::object: write_obj_line_(fs, "o", name); break;
    case obj_element::group: write_obj_line_(fs, "g", name); break;
    case obj_element::usemtl: write_obj_line_(fs, "usemtl", name); break;
    case obj_element::smoothing: write_obj_line_(fs, "s", name); break;
    case obj_element::mtllib: write_obj_line_(fs, "mtllib", name); break;
    case obj_element::objxlib: break;
  }
}
void write_mtl_element(
    FILE* fs, mtl_element element, const mtl_material& material) {
  switch (element) {
    case mtl_element::material: {
      static auto def    = mtl_material{};
      auto write_obj_opt = [](FILE* fs, const char* name, const auto& val,
                               const auto& def) {
        if (val == def) return;
        write_obj_line_(fs, name, val);
      };
      auto write_obj_txt = [](FILE* fs, const char* name,
                               const mtl_texture_info& info) {
        if (info.path.empty()) return;
        write_obj_line_(fs, name, info.path);
      };
      write_obj_line_(fs, "newmtl", material.name);
      write_obj_opt(fs, "  illum", material.illum, def.illum);
      write_obj_opt(fs, "  Ke", material.ke, def.ke);
      write_obj_opt(fs, "  Ka", material.ka, def.ka);
      write_obj_opt(fs, "  Kd", material.kd, vec3f{-1, -1, -1});
      write_obj_opt(fs, "  Ks", material.ks, vec3f{-1, -1, -1});
      write_obj_opt(fs, "  Kr", material.kr, def.kr);
      write_obj_opt(fs, "  Kt", material.kt, def.kt);
      write_obj_opt(fs, "  Ns", (int)material.ns, -1);
      write_obj_opt(fs, "  d", material.op, def.op);
      write_obj_opt(fs, "  Ni", material.ior, def.ior);
      write_obj_txt(fs, "  map_Ke", material.ke_map);
      write_obj_txt(fs, "  map_Ka", material.ka_map);
      write_obj_txt(fs, "  map_Kd", material.kd_map);
      write_obj_txt(fs, "  map_Ks", material.ks_map);
      write_obj_txt(fs, "  map_Kr", material.kr_map);
      write_obj_txt(fs, "  map_Kt", material.kt_map);
      write_obj_txt(fs, "  map_d", material.op_map);
      write_obj_txt(fs, "  map_Ni", material.ior_map);
      write_obj_txt(fs, "  map_bump", material.bump_map);
      write_obj_txt(fs, "  map_norm", material.norm_map);
      write_obj_txt(fs, "  map_disp", material.disp_map);
      write_obj_txt(fs, "  map_occ", material.occ_map);
      write_obj_opt(fs, "  Pr", material.pr, def.pr);
      write_obj_opt(fs, "  Pm", material.pm, def.pm);
      write_obj_opt(fs, "  Ps", material.ps, def.ps);
      write_obj_opt(fs, "  Pc", material.pc, def.pc);
      write_obj_opt(fs, "  Pcr", material.pcr, def.pcr);
      write_obj_txt(fs, "  Pr_map", material.pr_map);
      write_obj_txt(fs, "  Pm_map", material.pm_map);
      write_obj_txt(fs, "  Ps_map", material.ps_map);
      write_obj_txt(fs, "  Pc_map", material.pc_map);
      write_obj_txt(fs, "  Pcr_map", material.pcr_map);
      write_obj_opt(fs, "  Vt", material.vt, def.vt);
      write_obj_opt(fs, "  Ve", material.ve, def.ve);
      write_obj_opt(fs, "  Vs", material.vs, def.vs);
      write_obj_opt(fs, "  Vg", material.vg, def.vg);
      write_obj_opt(fs, "  Vr", material.vr, def.vr);
      write_obj_txt(fs, "  Vs_map", material.vs_map);
      write_obj_text(fs, "\n");
    } break;
  }
}
void write_objx_element(FILE* fs, objx_element element,
    const objx_camera& camera, const objx_environment& environment,
    const objx_instance& instance, const objx_procedural& procedural) {
  switch (element) {
    case objx_element::camera: {
      write_obj_line_(fs, "c", camera.name, (int)camera.ortho, camera.width,
          camera.height, camera.lens, camera.focus, camera.aperture,
          camera.frame);
    } break;
    case objx_element::environment: {
      write_obj_line_(fs, "e", environment.name, environment.ke,
          environment.ke_txt.path != "" ? environment.ke_txt.path : "\"\" "s,
          environment.frame);
    } break;
    case objx_element::instance: {
      write_obj_line_(fs, "i", instance.name, instance.object,
          instance.material, instance.frame);
    } break;
    case objx_element::procedural: {
      write_obj_line_(fs, "po", procedural.name, procedural.type,
          procedural.material, procedural.size, procedural.level,
          procedural.frame);
    } break;
  }
}

}  // namespace yocto
