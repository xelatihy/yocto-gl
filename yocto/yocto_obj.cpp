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
static inline void parse_obj_value(string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

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
  auto material = obj_material{};
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
      material = obj_material{};
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
      auto camera = obj_camera();
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
      auto environment = obj_environment();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.ke);
      parse_obj_value(line, environment.ke_txt.path);
      parse_obj_value(line, environment.frame);
      if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
      cb.environmnet(environment);
    } else if (cmd == "i") {
      auto instance = obj_instance();
      parse_obj_value(line, instance.name);
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.material);
      parse_obj_value(line, instance.frame);
      cb.instance(instance);
    } else if (cmd == "po") {
      auto procedural = obj_procedural();
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
    
void parse_obj_value(string_view& str, obj_value& value, obj_value_type type, int array_size = 3) {
    switch(type) {
        case obj_value_type::number: {
            auto value_ = 0.0f;
            parse_obj_value(str, value_);
            value = make_obj_value(value_);
        } break;
        case obj_value_type::string: {
            auto value_ = ""s;
            parse_obj_value(str, value_);
            value = make_obj_value(value_);
        } break;
        case obj_value_type::array: {
            if(array_size == 2) {
                auto value_ = zero2f;
                parse_obj_value(str, value_);
                value = make_obj_value(value_);
            } else
            if(array_size == 3) {
                auto value_ = zero3f;
                parse_obj_value(str, value_);
                value = make_obj_value(value_);
            } else
            if(array_size == 12) {
                auto value_ = identity3x4f;
                parse_obj_value(str, value_);
                value = make_obj_value(value_);
            } else {
                throw std::runtime_error("should not have gotten here");
            }
            } break;
        case obj_value_type::boolean: {
            auto value_ = 0;
            parse_obj_value(str, value_);
            value = make_obj_value((bool)value_);
        } break;
    }
}

static inline void parse_obj_value_or_empty(string_view& str, obj_value& value) {
    skip_obj_whitespace(str);
    if (str.empty()) {
        value = make_obj_value(""s);
    } else {
        parse_obj_value(str, value, obj_value_type::string);
    }
}

// Read obj
bool read_obj_command(FILE* fs, obj_command& command, obj_value& value,
    vector<obj_vertex>& vertices, obj_vertex& vert_size) {
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
      parse_obj_value(line, value, obj_value_type::array);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      parse_obj_value(line, value, obj_value_type::array);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
        parse_obj_value(line, value, obj_value_type::array, 2);
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
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
bool read_mtl_command(FILE* fs, mtl_command& command, obj_value& value,
    obj_texture_info& texture, bool fliptr) {
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
      command = mtl_command::material;
      parse_obj_value(line, value, obj_value_type::string);
    } else if (cmd == "illum") {
      command = mtl_command::illum;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Ke") {
      command = mtl_command::emission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Kd") {
      command = mtl_command::diffuse;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Ks") {
      command = mtl_command::specular;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Kt") {
      command = mtl_command::transmission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Tf") {
      command = mtl_command::transmission;
        auto color = vec3f{-1};
      value   = make_obj_value(color);
      parse_obj_value(line, value, obj_value_type::array);
        get_obj_value(value, color);
      if (color.y < 0) color = vec3f{color.x};
        if (fliptr) color = 1 - color;
        value   = make_obj_value(color);
    } else if (cmd == "Tr") {
      command = mtl_command::opacity;
      parse_obj_value(line, value, obj_value_type::number);
      if (fliptr) value.number = 1 - value.number;
    } else if (cmd == "Ns") {
      command = mtl_command::exponent;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "d") {
      command = mtl_command::opacity;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Ke") {
      command = mtl_command::emission_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Kd") {
      command = mtl_command::diffuse_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Ks") {
      command = mtl_command::specular_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Tr") {
      command = mtl_command::transmission_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      command = mtl_command::opacity_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_bump" || cmd == "bump") {
      command = mtl_command::bump_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_disp" || cmd == "disp") {
      command = mtl_command::displacement_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_norm" || cmd == "norm") {
      command = mtl_command::normal_map;
      parse_obj_value(line, texture);
    } else if (cmd == "Pm") {
      command = mtl_command::pbr_metallic;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pr") {
      command = mtl_command::pbr_roughness;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Ps") {
      command = mtl_command::pbr_sheen;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pc") {
      command = mtl_command::pbr_clearcoat;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pcr") {
      command = mtl_command::pbr_coatroughness;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Pm") {
      command = mtl_command::pbr_metallic_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pr") {
      command = mtl_command::pbr_roughness_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Ps") {
      command = mtl_command::pbr_sheen_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pc") {
      command = mtl_command::pbr_clearcoat_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pcr") {
      command = mtl_command::pbr_coatroughness_map;
      parse_obj_value(line, texture);
    } else if (cmd == "Vt") {
      command = mtl_command::vol_transmission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vp") {
      command = mtl_command::vol_meanfreepath;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Ve") {
      command = mtl_command::vol_emission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vs") {
      command = mtl_command::vol_scattering;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vg") {
      command = mtl_command::vol_anisotropy;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Vr") {
      command = mtl_command::vol_scale;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Vs") {
      command = mtl_command::vol_scattering_map;
      parse_obj_value(line, texture);
    } else {
      continue;
    }

    return true;
  }

  return false;
}

// Read objx
bool read_objx_command(FILE* fs, objx_command& command, obj_value& value,
    obj_texture_info& texture) {
  // read the file line by line
  char buffer[4096];
  auto pos = ftell(fs);
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

    // read values
    if (cmd == "newcam") {
      command = objx_command::camera;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newenv") {
      command = objx_command::environment;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newist") {
      command = objx_command::instance;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newproc") {
      command = objx_command::procedural;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "frame") {
      command = objx_command::frame;
      parse_obj_value(line, value, obj_value_type::array, 12);
      return true;
    } else if (cmd == "obj") {
      command = objx_command::object;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "mat") {
      command = objx_command::material;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "ortho") {
      command = objx_command::ortho;
      parse_obj_value(line, value, obj_value_type::boolean);
      return true;
    } else if (cmd == "width") {
      command = objx_command::width;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "height") {
      command = objx_command::height;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "lens") {
      command = objx_command::lens;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "aperture") {
      command = objx_command::aperture;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "focus") {
      command = objx_command::focus;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "Ke") {
      command = objx_command::emission;
      parse_obj_value(line, value, obj_value_type::array);
      return true;
    } else if (cmd == "map_Ke") {
      command = objx_command::emission_map;
      parse_obj_value(line, texture);
      return true;
    }
    // backward compatibility
    else if (cmd == "c") {
      auto oname    = value.string;
      auto name    =  obj_value{}, ortho    = obj_value{}, width    = obj_value{}, height   = obj_value{}, lens     = obj_value{}, aperture = obj_value{}, focus    = obj_value{}, frame    = obj_value{};
        parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, ortho, obj_value_type::boolean);
      parse_obj_value(line, width, obj_value_type::number);
      parse_obj_value(line, height, obj_value_type::number);
      parse_obj_value(line, lens, obj_value_type::number);
      parse_obj_value(line, focus, obj_value_type::number);
      parse_obj_value(line, aperture, obj_value_type::number);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::camera && oname != "") {
        command = objx_command::ortho;
        value   = ortho;
      } else if (command == objx_command::ortho) {
        command = objx_command::width;
        value   = width;
      } else if (command == objx_command::width) {
        command = objx_command::height;
        value   = height;
      } else if (command == objx_command::height) {
        command = objx_command::lens;
        value   = lens;
      } else if (command == objx_command::lens) {
        command = objx_command::focus;
        value   = focus;
      } else if (command == objx_command::focus) {
        command = objx_command::aperture;
        value   = aperture;
      } else if (command == objx_command::aperture) {
        command = objx_command::frame;
          value   = frame;
      } else {
        command = objx_command::camera;
          value   = name;
      }
      if (command != objx_command::frame) fseek(fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "e") {
        auto name = obj_value{}, frame = obj_value{}, emission = obj_value{}, emission_map = obj_value{};
        parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, emission, obj_value_type::array);
      parse_obj_value(line, emission_map, obj_value_type::string);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (emission_map.string == "\"\"") emission_map.string = "";
      if (command == objx_command::environment) {
        command = objx_command::emission;
          value   = emission;
      } else if (command == objx_command::emission) {
        command = objx_command::emission_map;
          get_obj_value(emission_map, texture.path);
      } else if (command == objx_command::emission_map) {
        command = objx_command::frame;
          value   = frame;
      } else {
        command = objx_command::environment;
          value   = name;
      }
      if (command != objx_command::frame) fseek(fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "i") {
        auto name = obj_value{}, frame = obj_value{}, object = obj_value{}, material = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, object, obj_value_type::string);
      parse_obj_value(line, material, obj_value_type::string);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::instance) {
        command = objx_command::object;
          value   = object;
      } else if (command == objx_command::object) {
        command = objx_command::material;
          value   = material;
      } else if (command == objx_command::material) {
        command = objx_command::frame;
          value   = frame;
      } else {
        command = objx_command::instance;
          value   = name;
      }
      if (command != objx_command::frame) fseek(fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "po") {
        auto name = obj_value{}, frame = obj_value{}, type = obj_value{}, material = obj_value{}, size = obj_value{}, level = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, type, obj_value_type::string);
      parse_obj_value(line, material, obj_value_type::string);
      parse_obj_value(line, size, obj_value_type::number);
      parse_obj_value(line, level, obj_value_type::number);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::procedural) {
        command = objx_command::object;
        value    = type;
      } else if (command == objx_command::object) {
        command = objx_command::material;
        value    = material;
      } else if (command == objx_command::material) {
        command = objx_command::frame;
          value = frame;
      } else {
        command = objx_command::procedural;
          value = name;
      }
      if (command != objx_command::frame) fseek(fs, pos, SEEK_SET);
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

void write_obj_command(FILE* fs, obj_command command, const obj_value& value_,
    const vector<obj_vertex>& vertices) {
    auto& name = value_.string;
    auto& value = value_.array_;
  switch (command) {
    case obj_command::vertex:
      checked_fprintf(fs, "v %g %g %g\n", value[0], value[1], value[2]);
      break;
    case obj_command::normal:
      checked_fprintf(fs, "vn %g  %g %g\n", value[0], value[1], value[2]);
      break;
    case obj_command::texcoord:
      checked_fprintf(fs, "vt %g %g\n", value[0], value[1]);
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

void write_mtl_command(FILE* fs, mtl_command command, const obj_value& value_,
    const obj_texture_info& texture) {
    auto& name = value_.string;
    auto value = value_.number;
    auto& color = value_.array_;
  switch (command) {
    case mtl_command::material:
      checked_fprintf(fs, "\nnewmtl %s\n", name.c_str());
      break;
    case mtl_command::illum:
      checked_fprintf(fs, "  illum %d\n", (int)value);
      break;
    case mtl_command::emission:
      checked_fprintf(fs, "  Ke %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::ambient:
      checked_fprintf(fs, "  Ka %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::diffuse:
      checked_fprintf(fs, "  Kd %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::specular:
      checked_fprintf(fs, "  Ks %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::reflection:
      checked_fprintf(fs, "  Kr %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::transmission:
      checked_fprintf(fs, "  Kt %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::exponent:
      checked_fprintf(fs, "  Ns %d\n", (int)value);
      break;
    case mtl_command::opacity: checked_fprintf(fs, "  d %g\n", value); break;
    case mtl_command::ior: checked_fprintf(fs, "  Ni %g\n", value); break;
    case mtl_command::emission_map:
      checked_fprintf(fs, "  map_Ke %s\n", texture.path.c_str());
      break;
    case mtl_command::ambient_map:
      checked_fprintf(fs, "  map_Ka %s\n", texture.path.c_str());
      break;
    case mtl_command::diffuse_map:
      checked_fprintf(fs, "  map_Kd %s\n", texture.path.c_str());
      break;
    case mtl_command::specular_map:
      checked_fprintf(fs, "  map_Ks %s\n", texture.path.c_str());
      break;
    case mtl_command::reflection_map:
      checked_fprintf(fs, "  map_Kr %s\n", texture.path.c_str());
      break;
    case mtl_command::transmission_map:
      checked_fprintf(fs, "  map_Kt %s\n", texture.path.c_str());
      break;
    case mtl_command::opacity_map:
      checked_fprintf(fs, "  map_d %s\n", texture.path.c_str());
      break;
    case mtl_command::exponent_map:
      checked_fprintf(fs, "  map_Ni %s\n", texture.path.c_str());
      break;
    case mtl_command::bump_map:
      checked_fprintf(fs, "  map_bump %s\n", texture.path.c_str());
      break;
    case mtl_command::normal_map:
      checked_fprintf(fs, "  map_norm %s\n", texture.path.c_str());
      break;
    case mtl_command::displacement_map:
      checked_fprintf(fs, "  map_disp %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_roughness:
      checked_fprintf(fs, "  Pr %g\n", value);
      break;
    case mtl_command::pbr_metallic:
      checked_fprintf(fs, "  Pm %g\n", value);
      break;
    case mtl_command::pbr_sheen: checked_fprintf(fs, "  Ps %g\n", value); break;
    case mtl_command::pbr_clearcoat:
      checked_fprintf(fs, "  Pc %g\n", value);
      break;
    case mtl_command::pbr_coatroughness:
      checked_fprintf(fs, "  Pcr %g\n", value);
      break;
    case mtl_command::pbr_roughness_map:
      checked_fprintf(fs, "  Pr_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_metallic_map:
      checked_fprintf(fs, "  Pm_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_sheen_map:
      checked_fprintf(fs, "  Ps_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_clearcoat_map:
      checked_fprintf(fs, "  Pc_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_coatroughness_map:
      checked_fprintf(fs, "  Pcr_map %s\n", texture.path.c_str());
      break;
    case mtl_command::vol_transmission:
      checked_fprintf(fs, "  Vt %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_meanfreepath:
      checked_fprintf(fs, "  Vp %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_emission:
      checked_fprintf(fs, "  Ve %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_scattering:
      checked_fprintf(fs, "  Vs %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_anisotropy:
      checked_fprintf(fs, "  Vg %g\n", value);
      break;
    case mtl_command::vol_scale: checked_fprintf(fs, "  Vr %g\n", value); break;
    case mtl_command::vol_scattering_map:
      checked_fprintf(fs, "  Vs_map %s\n", texture.path.c_str());
  }
}

void write_objx_command(FILE* fs, objx_command command, const obj_value& value_,
    const obj_texture_info& texture) {
    auto& name = value_.string;
    auto value = value_.number;
    auto& color = value_.array_;
    auto& frame = value_.array_;
  switch (command) {
    case objx_command::camera:
      checked_fprintf(fs, "\nnewcam %s\n", name.c_str());
      break;
    case objx_command::environment:
      checked_fprintf(fs, "\nnewenv %s\n", name.c_str());
      break;
    case objx_command::instance:
      checked_fprintf(fs, "\nnewist %s\n", name.c_str());
      break;
    case objx_command::procedural:
      checked_fprintf(fs, "\nnewproc %s\n", name.c_str());
      break;
    case objx_command::frame:
      checked_fprintf(fs, "  frame %g %g %g %g %g %g %g %g %g %g %g %g\n",
          frame[0], frame[1], frame[2], frame[3], frame[4], frame[5],
          frame[6], frame[7], frame[8], frame[9], frame[10], frame[11]);
      break;
    case objx_command::object:
      checked_fprintf(fs, "  obj %s\n", name.c_str());
      break;
    case objx_command::material:
      checked_fprintf(fs, "  mat %s\n", name.c_str());
      break;
    case objx_command::ortho: checked_fprintf(fs, "  ortho %g\n", value); break;
    case objx_command::width: checked_fprintf(fs, "  width %g\n", value); break;
    case objx_command::height:
      checked_fprintf(fs, "  height %g\n", value);
      break;
    case objx_command::lens: checked_fprintf(fs, "  lens %g\n", value); break;
    case objx_command::aperture:
      checked_fprintf(fs, "  aperture %g\n", value);
      break;
    case objx_command::focus: checked_fprintf(fs, "  focus %g\n", value); break;
    case objx_command::emission:
      checked_fprintf(fs, "  Ke %g %g %g\n", color[0], color[1], color[2]);
      break;
    case objx_command::emission_map:
      checked_fprintf(fs, "  map_Ke %s\n", texture.path.c_str());
      break;
  }
}

}  // namespace yocto
