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

    // get command
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

    // get command
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

    // get command
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
bool read_obj_command(FILE* fs, obj_command& command, vec3f& value,
    string& name, vector<obj_vertex>& vertices, obj_vertex& vert_size) {
  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      parse_obj_value(line, value);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      parse_obj_value(line, value);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
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
      if (cmd == "f") command = obj_command::face;
      if (cmd == "l") command = obj_command::line;
      if (cmd == "p") command = obj_command::point;
      return true;
    } else if (cmd == "o") {
      command = obj_command::object;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      parse_obj_value(line, name);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
bool read_mtl_command(
    FILE* fs, mtl_command& command, mtl_material& material, bool fliptr) {
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

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      if (found) {
        fseek(fs, fpos, SEEK_SET);
        command = mtl_command::material;
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
    command = mtl_command::material;
    return true;
  } else {
    return false;
  }
}

// Read objx
bool read_objx_command(FILE* fs, objx_command& command, objx_camera& camera,
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

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "c") {
      command = objx_command::camera;
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
      command     = objx_command::environment;
      environment = objx_environment();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.ke);
      parse_obj_value(line, environment.ke_txt.path);
      parse_obj_value(line, environment.frame);
      if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
      return true;
    } else if (cmd == "i") {
      command  = objx_command::instance;
      instance = objx_instance();
      parse_obj_value(line, instance.name);
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.material);
      parse_obj_value(line, instance.frame);
      return true;
    } else if (cmd == "po") {
      command    = objx_command::procedural;
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

static inline void checked_fprintf(FILE* fs, const char* fmt, ...) {
  va_list args1;
  va_start(args1, fmt);
  if (vfprintf(fs, fmt, args1) < 0)
    throw std::runtime_error("cannot write to file");
  va_end(args1);
}

// Write obj elements
void write_obj_comment(FILE* fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    checked_fprintf(fs, "# %s\n", line.c_str());
  }
  checked_fprintf(fs, "\n");
}

void write_obj_command(FILE* fs, obj_command command, const vec3f& value,
    const string& name, const vector<obj_vertex>& vertices) {
  switch (command) {
    case obj_command::vertex:
      checked_fprintf(fs, "v %g %g %g\n", value.x, value.y, value.z);
      break;
    case obj_command::normal:
      checked_fprintf(fs, "vn %g  %g %g\n", value.x, value.y, value.z);
      break;
    case obj_command::texcoord:
      checked_fprintf(fs, "vt %g %g\n", value.x, value.y);
      break;
    case obj_command::face:
    case obj_command::line:
    case obj_command::point:
      if (command == obj_command::face) checked_fprintf(fs, "f ");
      if (command == obj_command::line) checked_fprintf(fs, "l ");
      if (command == obj_command::point) checked_fprintf(fs, "p ");
      for (auto& vert : vertices) {
        checked_fprintf(fs, " ");
        checked_fprintf(fs, "%d", vert.position);
        if (vert.texcoord) {
          checked_fprintf(fs, "/%d", vert.texcoord);
          if (vert.normal) {
            checked_fprintf(fs, "/%d", vert.normal);
          }
        } else if (vert.normal) {
          checked_fprintf(fs, "//%d", vert.normal);
        }
      }
      checked_fprintf(fs, "\n");
      break;
    case obj_command::object:
      checked_fprintf(fs, "o %s\n", name.c_str());
      break;
    case obj_command::group: checked_fprintf(fs, "g %s\n", name.c_str()); break;
    case obj_command::usemtl:
      checked_fprintf(fs, "usemtl %s\n", name.c_str());
      break;
    case obj_command::smoothing:
      checked_fprintf(fs, "s %s\n", name.c_str());
      break;
    case obj_command::mtllib:
      checked_fprintf(fs, "mtllib %s\n", name.c_str());
      break;
    case obj_command::objxlib: break;
  }
}
void write_mtl_command(
    FILE* fs, mtl_command command, const mtl_material& material) {
  switch (command) {
    case mtl_command::material: {
      static auto def = mtl_material{};
      checked_fprintf(fs, "newmtl %s\n", material.name.c_str());
      checked_fprintf(fs, "  illum %d\n", material.illum);
      if (material.ke != def.ke)
        checked_fprintf(
            fs, "  Ke %g %g %g\n", material.ke.x, material.ke.y, material.ke.z);
      if (material.ka != def.ka)
        checked_fprintf(
            fs, "  Ka %g %g %g\n", material.ka.x, material.ka.y, material.ka.z);
      checked_fprintf(
          fs, "  Kd %g %g %g\n", material.kd.x, material.kd.y, material.kd.z);
      checked_fprintf(
          fs, "  Ks %g %g %g\n", material.ks.x, material.ks.y, material.ks.z);
      if (material.kr != def.kr)
        checked_fprintf(
            fs, "  Kr %g %g %g\n", material.kr.x, material.kr.y, material.kr.z);
      if (material.kt != def.kt)
        checked_fprintf(
            fs, "  Kt %g %g %g\n", material.kt.x, material.kt.y, material.kt.z);
      checked_fprintf(fs, "  Ns %d\n", (int)material.ns, -1);
      if (material.op != def.op) checked_fprintf(fs, "  d %g\n", material.op);
      if (material.ior != def.ior)
        checked_fprintf(fs, "  Ni %g\n", material.ior);
      if (!material.ke_map.path.empty())
        checked_fprintf(fs, "  map_Ke %s\n", material.ke_map.path.c_str());
      if (!material.ka_map.path.empty())
        checked_fprintf(fs, "  map_Ka %s\n", material.ka_map.path.c_str());
      if (!material.kd_map.path.empty())
        checked_fprintf(fs, "  map_Kd %s\n", material.kd_map.path.c_str());
      if (!material.ks_map.path.empty())
        checked_fprintf(fs, "  map_Ks %s\n", material.ks_map.path.c_str());
      if (!material.kr_map.path.empty())
        checked_fprintf(fs, "  map_Kr %s\n", material.kr_map.path.c_str());
      if (!material.kt_map.path.empty())
        checked_fprintf(fs, "  map_Kt %s\n", material.kt_map.path.c_str());
      if (!material.op_map.path.empty())
        checked_fprintf(fs, "  map_d %s\n", material.op_map.path.c_str());
      if (!material.ior_map.path.empty())
        checked_fprintf(fs, "  map_Ni %s\n", material.ior_map.path.c_str());
      if (!material.bump_map.path.empty())
        checked_fprintf(fs, "  map_bump %s\n", material.bump_map.path.c_str());
      if (!material.norm_map.path.empty())
        checked_fprintf(fs, "  map_norm %s\n", material.norm_map.path.c_str());
      if (!material.disp_map.path.empty())
        checked_fprintf(fs, "  map_disp %s\n", material.disp_map.path.c_str());
      if (!material.occ_map.path.empty())
        checked_fprintf(fs, "  map_occ %s\n", material.occ_map.path.c_str());
      if (material.pr != def.pr) checked_fprintf(fs, "  Pr %g\n", material.pr);
      if (material.pm != def.pm) checked_fprintf(fs, "  Pm %g\n", material.pm);
      if (material.ps != def.ps) checked_fprintf(fs, "  Ps %g\n", material.ps);
      if (material.pc != def.pc) checked_fprintf(fs, "  Pc %g\n", material.pc);
      if (material.pcr != def.pcr)
        checked_fprintf(fs, "  Pcr %g\n", material.pcr);
      if (!material.pr_map.path.empty())
        checked_fprintf(fs, "  Pr_map %s\n", material.pr_map.path.c_str());
      if (!material.pm_map.path.empty())
        checked_fprintf(fs, "  Pm_map %s\n", material.pm_map.path.c_str());
      if (!material.ps_map.path.empty())
        checked_fprintf(fs, "  Ps_map %s\n", material.ps_map.path.c_str());
      if (!material.pc_map.path.empty())
        checked_fprintf(fs, "  Pc_map %s\n", material.pc_map.path.c_str());
      if (!material.pcr_map.path.empty())
        checked_fprintf(fs, "  Pcr_map %s\n", material.pcr_map.path.c_str());
      if (material.vt != def.vt)
        checked_fprintf(
            fs, "  Vt %g %g %g\n", material.vt.x, material.vt.y, material.vt.z);
      if (material.ve != def.ve)
        checked_fprintf(
            fs, "  Ve %g %g %g\n", material.ve.x, material.ve.y, material.ve.z);
      if (material.vs != def.vs)
        checked_fprintf(
            fs, "  Vs %g %g %g\n", material.vs.x, material.vs.y, material.vs.z);
      if (material.vg != def.vg) checked_fprintf(fs, "  Vg %g\n", material.vg);
      if (material.vr != def.vr) checked_fprintf(fs, "  Vr %g\n", material.vr);
      if (!material.vs_map.path.empty())
        checked_fprintf(fs, "  Vs_map %s\n", material.vs_map.path.c_str());
      checked_fprintf(fs, "\n");
    } break;
  }
}
void write_objx_command(FILE* fs, objx_command command,
    const objx_camera& camera, const objx_environment& environment,
    const objx_instance& instance, const objx_procedural& procedural) {
  switch (command) {
    case objx_command::camera: {
      checked_fprintf(fs,
          "c %s %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
          camera.name.c_str(), (int)camera.ortho, camera.width, camera.height,
          camera.lens, camera.focus, camera.aperture, camera.frame.x.x,
          camera.frame.x.y, camera.frame.x.z, camera.frame.y.x,
          camera.frame.y.y, camera.frame.y.z, camera.frame.z.x,
          camera.frame.z.y, camera.frame.z.z, camera.frame.o.x,
          camera.frame.o.y, camera.frame.o.z);
    } break;
    case objx_command::environment: {
      checked_fprintf(fs,
          "e %s %g %g %g %s %g %g %g %g %g %g %g %g %g %g %g %g\n",
          environment.name.c_str(), environment.ke.x, environment.ke.y,
          environment.ke.z,
          environment.ke_txt.path != "" ? environment.ke_txt.path.c_str()
                                        : "\"\" ",
          environment.frame.x.x, environment.frame.x.y, environment.frame.x.z,
          environment.frame.y.x, environment.frame.y.y, environment.frame.y.z,
          environment.frame.z.x, environment.frame.z.y, environment.frame.z.z,
          environment.frame.o.x, environment.frame.o.y, environment.frame.o.z);
    } break;
    case objx_command::instance: {
      checked_fprintf(fs, "i %s %s %s %g %g %g %g %g %g %g %g %g %g %g %g\n",
          instance.name.c_str(), instance.object.c_str(),
          instance.material.c_str(), instance.frame.x.x, instance.frame.x.y,
          instance.frame.x.z, instance.frame.y.x, instance.frame.y.y,
          instance.frame.y.z, instance.frame.z.x, instance.frame.z.y,
          instance.frame.z.z, instance.frame.o.x, instance.frame.o.y,
          instance.frame.o.z);
    } break;
    case objx_command::procedural: {
      checked_fprintf(fs,
          "po %s %s %s %f %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
          procedural.name.c_str(), procedural.type.c_str(),
          procedural.material.c_str(), procedural.size, procedural.level,
          procedural.frame.x.x, procedural.frame.x.y, procedural.frame.x.z,
          procedural.frame.y.x, procedural.frame.y.y, procedural.frame.y.z,
          procedural.frame.z.x, procedural.frame.z.y, procedural.frame.z.z,
          procedural.frame.o.x, procedural.frame.o.y, procedural.frame.o.z);
    } break;
  }
}

}  // namespace yocto
