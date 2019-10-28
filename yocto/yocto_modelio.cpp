//
// Implementation for Yocto/ModelIO.
//

//
// TODO: remove obj procedurals
// TODO: simplify objx parsing using single line commands
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

#include "yocto_modelio.h"

#include "yocto_commonio.h"
#include "yocto_image.h"

#include <algorithm>
#include <cinttypes>
#include <climits>
#include <cstdarg>
#include <limits>
#include <string_view>

#define CGLTF_IMPLEMENTATION
#include "ext/cgltf.h"

// -----------------------------------------------------------------------------
// FILE AND PROPERTY HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
file_wrapper::file_wrapper(file_wrapper&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
file_wrapper::~file_wrapper() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
void open_file(file_wrapper& fs, const string& filename, const string& mode) {
  close_file(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
file_wrapper open_file(const string& filename, const string& mode) {
  auto fs = file_wrapper{};
  open_file(fs, filename, mode);
  return fs;
}
void close_file(file_wrapper& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;
using namespace std::literals::string_view_literals;

template <typename T>
static inline T swap_endian(T value) {
  // https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
  union {
    T             value;
    unsigned char bytes[sizeof(T)];
  } source, dest;
  source.value = value;
  for (auto k = (size_t)0; k < sizeof(T); k++)
    dest.bytes[k] = source.bytes[sizeof(T) - k - 1];
  return dest.value;
}

// Read a line
static inline bool read_line(file_wrapper& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

static inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
static inline bool is_digit(char c) { return c >= '0' && c <= '9'; }
static inline bool is_alpha(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

static inline void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}
static inline void trim_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
  while (!str.empty() && is_space(str.back())) str.remove_suffix(1);
}

static inline bool is_whitespace(string_view str) {
  while (!str.empty()) {
    if (!is_space(str.front())) return false;
    str.remove_prefix(1);
  }
  return true;
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

static inline vector<vec2f> flip_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Parse values from a string
static inline void parse_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
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
static inline void parse_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_value(str, valuev);
  value = string{valuev};
}
static inline void parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_value(str, valuei);
  value = (bool)valuei;
}
static inline void parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
static inline void parse_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_value(str, values[i]);
}

static inline void parse_value(string_view& str, vec2f& value) {
  parse_value(str, &value.x, 2);
}
static inline void parse_value(string_view& str, vec3f& value) {
  parse_value(str, &value.x, 3);
}
static inline void parse_value(string_view& str, frame3f& value) {
  parse_value(str, &value.x.x, 12);
}
#ifdef __APPLE__
static inline void parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
#endif

// Parse values from a string
template <typename T>
static inline void parse_value_or_empty(string_view& str, T& value) {
  skip_whitespace(str);
  if (str.empty()) {
    value = T{};
  } else {
    parse_value(str, value);
  }
}

// Formats values to string
static inline void format_value(string& str, const string& value) {
  str += value;
}
static inline void format_value(string& str, const char* value) {
  str += value;
}
static inline void format_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
static inline void format_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
static inline void format_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static inline void format_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static inline void format_value(string& str, const vec2f& value) {
  char buf[256];
  sprintf(buf, "%g %g", value.x, value.y);
  str += buf;
}
static inline void format_value(string& str, const vec3f& value) {
  char buf[256];
  sprintf(buf, "%g %g %g", value.x, value.y, value.z);
  str += buf;
}
#if 0
static inline void format_value(string& str, const vec2i& value) {
  char buf[256];
  sprintf(buf, "%d %d", value.x, value.y);
  str += buf;
}
static inline void format_value(string& str, const vec3i& value) {
  char buf[256];
  sprintf(buf, "%d %d %d", value.x, value.y, value.z);
  str += buf;
}
#endif
static inline void format_value(string& str, const frame3f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g", value.x.x, value.x.y,
      value.x.z, value.y.x, value.y.y, value.y.z, value.z.x, value.z.y,
      value.z.z, value.o.x, value.o.y, value.o.z);
  str += buf;
}
static inline void format_value(string& str, const mat4f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", value.x.x,
      value.x.y, value.x.z, value.x.w, value.y.x, value.y.y, value.y.z,
      value.y.w, value.z.x, value.z.y, value.z.z, value.z.w, value.w.x,
      value.w.y, value.w.z, value.w.w);
  str += buf;
}

// Foramt to file
static inline void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static inline void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
static inline void format_values(
    file_wrapper& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
static inline void format_value(file_wrapper& fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

static inline void write_text(file_wrapper& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
static inline void write_text(file_wrapper& fs, const char* value) {
  if (fputs(value, fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}

template <typename T>
static inline void write_value(file_wrapper& fs, const T& value) {
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}
template <typename T>
static inline void write_value(
    file_wrapper& fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}

template <typename T>
static inline void read_value(file_wrapper& fs, T& value) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
}
template <typename T>
static inline void read_value(file_wrapper& fs, T& value, bool big_endian) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
  if (big_endian) value = swap_endian(value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE GLTF LOADER
// -----------------------------------------------------------------------------
namespace yocto {

void update_transforms(
    gltf_model& scene, gltf_node& node, const frame3f& parent = identity3x4f) {
  auto frame = parent * node.local * translation_frame(node.translation) *
               rotation_frame(node.rotation) * scaling_frame(node.scale);
  for (auto child : node.children)
    update_transforms(scene, scene.nodes[child], frame);
}

// convert gltf to scene
void load_gltf(const string& filename, gltf_model& scene) {
  // load gltf
  auto params = cgltf_options{};
  memset(&params, 0, sizeof(params));
  auto data   = (cgltf_data*)nullptr;
  auto result = cgltf_parse_file(&params, filename.c_str(), &data);
  if (result != cgltf_result_success) {
    throw std::runtime_error("could not load " + filename);
  }
  auto gltf = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      data, cgltf_free};
  auto dirname = get_dirname(filename);
  if (dirname != "") dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success) {
    throw std::runtime_error("could not load gltf buffers " + filename);
  }

  // convert textures
  auto _startswith = [](string_view str, string_view substr) {
    if (str.size() < substr.size()) return false;
    return str.substr(0, substr.size()) == substr;
  };
  auto imap = unordered_map<cgltf_image*, int>{};
  for (auto tid = 0; tid < gltf->images_count; tid++) {
    auto gimg        = &gltf->images[tid];
    auto texture     = gltf_texture{};
    texture.name     = gimg->name ? gimg->name : "";
    texture.filename = (_startswith(gimg->uri, "data:"))
                           ? string("[glTF-static inline].png")
                           : gimg->uri;
    scene.textures.push_back(texture);
    imap[gimg] = tid;
  }

  // add a texture
  auto get_texture = [&imap](const cgltf_texture_view& ginfo) {
    if (!ginfo.texture || !ginfo.texture->image) return -1;
    auto gtxt = ginfo.texture;
    return imap.at(gtxt->image);
  };

  // convert materials
  auto mmap = unordered_map<cgltf_material*, int>{{nullptr, -1}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto gmat             = &gltf->materials[mid];
    mmap[gmat]            = mid;
    auto& material        = scene.materials.emplace_back();
    material.name         = gmat->name ? gmat->name : "";
    material.emission     = {gmat->emissive_factor[0], gmat->emissive_factor[1],
        gmat->emissive_factor[2]};
    material.emission_tex = get_texture(gmat->emissive_texture);
    if (gmat->has_pbr_specular_glossiness) {
      material.has_specgloss = true;
      auto gsg               = &gmat->pbr_specular_glossiness;
      material.sg_diffuse    = vec4f{gsg->diffuse_factor[0],
          gsg->diffuse_factor[1], gsg->diffuse_factor[2],
          gsg->diffuse_factor[3]};
      material.sg_specular = {gsg->specular_factor[0], gsg->specular_factor[1],
          gsg->specular_factor[2]};
      material.sg_glossiness   = gsg->glossiness_factor;
      material.sg_diffuse_tex  = get_texture(gsg->diffuse_texture);
      material.sg_specular_tex = get_texture(gsg->specular_glossiness_texture);
    } else if (gmat->has_pbr_metallic_roughness) {
      material.has_metalrough  = true;
      auto gmr                 = &gmat->pbr_metallic_roughness;
      material.mr_base         = vec4f{gmr->base_color_factor[0],
          gmr->base_color_factor[1], gmr->base_color_factor[2],
          gmr->base_color_factor[3]};
      material.mr_metallic     = gmr->metallic_factor;
      material.mr_roughness    = gmr->roughness_factor;
      material.mr_base_tex     = get_texture(gmr->base_color_texture);
      material.mr_metallic_tex = get_texture(gmr->metallic_roughness_texture);
    }
    material.normal_tex = get_texture(gmat->normal_texture);
  }

  // get values from accessors
  auto accessor_values =
      [](const cgltf_accessor* gacc,
          bool normalize = false) -> vector<std::array<double, 4>> {
    auto gview       = gacc->buffer_view;
    auto data        = (byte*)gview->buffer->data;
    auto offset      = gacc->offset + gview->offset;
    auto stride      = gview->stride;
    auto compTypeNum = gacc->component_type;
    auto count       = gacc->count;
    auto type        = gacc->type;
    auto ncomp       = 0;
    if (type == cgltf_type_scalar) ncomp = 1;
    if (type == cgltf_type_vec2) ncomp = 2;
    if (type == cgltf_type_vec3) ncomp = 3;
    if (type == cgltf_type_vec4) ncomp = 4;
    auto compSize = 1;
    if (compTypeNum == cgltf_component_type_r_16 ||
        compTypeNum == cgltf_component_type_r_16u) {
      compSize = 2;
    }
    if (compTypeNum == cgltf_component_type_r_32u ||
        compTypeNum == cgltf_component_type_r_32f) {
      compSize = 4;
    }
    if (!stride) stride = compSize * ncomp;
    auto vals = vector<std::array<double, 4>>(count, {{0.0, 0.0, 0.0, 1.0}});
    for (auto i = 0; i < count; i++) {
      auto d = data + offset + i * stride;
      for (auto c = 0; c < ncomp; c++) {
        if (compTypeNum == cgltf_component_type_r_8) {  // char
          vals[i][c] = (double)(*(char*)d);
          if (normalize) vals[i][c] /= SCHAR_MAX;
        } else if (compTypeNum == cgltf_component_type_r_8u) {  // byte
          vals[i][c] = (double)(*(byte*)d);
          if (normalize) vals[i][c] /= UCHAR_MAX;
        } else if (compTypeNum == cgltf_component_type_r_16) {  // short
          vals[i][c] = (double)(*(short*)d);
          if (normalize) vals[i][c] /= SHRT_MAX;
        } else if (compTypeNum ==
                   cgltf_component_type_r_16u) {  // unsigned short
          vals[i][c] = (double)(*(unsigned short*)d);
          if (normalize) vals[i][c] /= USHRT_MAX;
        } else if (compTypeNum == cgltf_component_type_r_32u) {  // unsigned int
          vals[i][c] = (double)(*(unsigned int*)d);
          if (normalize) vals[i][c] /= UINT_MAX;
        } else if (compTypeNum == cgltf_component_type_r_32f) {  // float
          vals[i][c] = (*(float*)d);
        }
        d += compSize;
      }
    }
    return vals;
  };

  // convert meshes
  auto smap = unordered_map<cgltf_mesh*, int>{{nullptr, -1}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh  = &gltf->meshes[mid];
    smap[gmesh] = mid;
    auto& mesh  = scene.meshes.emplace_back();
    mesh.name   = gmesh->name ? gmesh->name : "";
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (!gprim->attributes_count) continue;
      auto& shape = mesh.primitives.emplace_back();
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name ? gattr->name : "");
        auto gacc     = gattr->data;
        auto vals     = accessor_values(gacc);
        if (semantic == "POSITION") {
          shape.positions.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.positions.push_back(
                {(float)vals[i][0], (float)vals[i][1], (float)vals[i][2]});
        } else if (semantic == "NORMAL") {
          shape.normals.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.normals.push_back(
                {(float)vals[i][0], (float)vals[i][1], (float)vals[i][2]});
        } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
          shape.texcoords.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.texcoords.push_back({(float)vals[i][0], (float)vals[i][1]});
        } else if (semantic == "COLOR" || semantic == "COLOR_0") {
          shape.colors.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.colors.push_back({(float)vals[i][0], (float)vals[i][1],
                (float)vals[i][2], (float)vals[i][3]});
        } else if (semantic == "TANGENT") {
          shape.tangents.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.tangents.push_back({(float)vals[i][0], (float)vals[i][1],
                (float)vals[i][2], (float)vals[i][3]});
          for (auto& t : shape.tangents) t.w = -t.w;
        } else if (semantic == "RADIUS") {
          shape.radius.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.radius.push_back((float)vals[i][0]);
        } else {
          // ignore
        }
      }
      // indices
      if (!gprim->indices) {
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.reserve(shape.positions.size() / 3);
          for (auto i = 0; i < shape.positions.size() / 3; i++)
            shape.triangles.push_back({i * 3 + 0, i * 3 + 1, i * 3 + 2});
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.reserve(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles.push_back({0, i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.reserve(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles.push_back({i - 2, i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.reserve(shape.positions.size() / 2);
          for (auto i = 0; i < shape.positions.size() / 2; i++)
            shape.lines.push_back({i * 2 + 0, i * 2 + 1});
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.reserve(shape.positions.size());
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines.push_back({i - 1, i});
          shape.lines.back() = {(int)shape.positions.size() - 1, 0};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.reserve(shape.positions.size() - 1);
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines.push_back({i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          throw std::runtime_error("points not supported");
        } else {
          throw std::runtime_error("unknown primitive type");
        }
      } else {
        auto indices = accessor_values(gprim->indices);
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.reserve(indices.size() / 3);
          for (auto i = 0; i < indices.size() / 3; i++)
            shape.triangles.push_back({(int)indices[i * 3 + 0][0],
                (int)indices[i * 3 + 1][0], (int)indices[i * 3 + 2][0]});
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.reserve(indices.size() - 2);
          for (auto i = 2; i < indices.size(); i++)
            shape.triangles.push_back({(int)indices[0][0],
                (int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.reserve(indices.size() - 2);
          for (auto i = 2; i < indices.size(); i++)
            shape.triangles.push_back({(int)indices[i - 2][0],
                (int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.reserve(indices.size() / 2);
          for (auto i = 0; i < indices.size() / 2; i++)
            shape.lines.push_back(
                {(int)indices[i * 2 + 0][0], (int)indices[i * 2 + 1][0]});
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.reserve(indices.size());
          for (auto i = 1; i < indices.size(); i++)
            shape.lines.push_back({(int)indices[i - 1][0], (int)indices[i][0]});
          shape.lines.back() = {
              (int)indices[indices.size() - 1][0], (int)indices[0][0]};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.reserve(indices.size() - 1);
          for (auto i = 1; i < indices.size(); i++)
            shape.lines.push_back({(int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_points) {
          throw std::runtime_error("points not supported");
        } else {
          throw std::runtime_error("unknown primitive type");
        }
      }
    }
  }

  // convert cameras
  auto cmap = unordered_map<cgltf_camera*, int>{{nullptr, -1}};
  for (auto cid = 0; cid < gltf->cameras_count; cid++) {
    auto gcam    = &gltf->cameras[cid];
    cmap[gcam]   = cid;
    auto& camera = scene.cameras.emplace_back();
    camera.name  = gcam->name ? gcam->name : "";
    camera.ortho = gcam->type == cgltf_camera_type_orthographic;
    if (camera.ortho) {
      // throw std::runtime_error("orthographic not supported well");
      auto ortho    = &gcam->orthographic;
      camera.yfov   = ortho->ymag;
      camera.aspect = ortho->xmag / ortho->ymag;
    } else {
      auto persp    = &gcam->perspective;
      camera.yfov   = persp->yfov;
      camera.aspect = persp->aspect_ratio;
    }
    scene.cameras.push_back(camera);
    cmap[gcam] = (int)scene.cameras.size() - 1;
  }

  // convert nodes
  auto nmap = unordered_map<cgltf_node*, int>{{nullptr, -1}};
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde  = &gltf->nodes[nid];
    nmap[gnde] = nid;
    auto& node = scene.nodes.emplace_back();
    node.name  = gnde->name ? gnde->name : "";
    if (gnde->camera) node.camera = cmap.at(gnde->camera);
    if (gnde->mesh) node.mesh = smap.at(gnde->mesh);
    if (gnde->has_translation) {
      node.translation = {
          gnde->translation[0], gnde->translation[1], gnde->translation[2]};
    }
    if (gnde->has_rotation) {
      node.rotation = {gnde->rotation[0], gnde->rotation[1], gnde->rotation[2],
          gnde->rotation[3]};
    }
    if (gnde->has_scale) {
      node.scale = {gnde->scale[0], gnde->scale[1], gnde->scale[2]};
    }
    if (gnde->has_matrix) {
      auto m     = gnde->matrix;
      node.local = frame3f(
          mat4f{{m[0], m[1], m[2], m[3]}, {m[4], m[5], m[6], m[7]},
              {m[8], m[9], m[10], m[11]}, {m[12], m[13], m[14], m[15]}});
    }
  }

  // set up parent pointers
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde = &gltf->nodes[nid];
    if (!gnde->children_count) continue;
    for (auto cid = 0; cid < gnde->children_count; cid++) {
      scene.nodes[nid].children.push_back(nmap.at(gnde->children[cid]));
      scene.nodes[nmap.at(gnde->children[cid])].parent = nid;
    }
  }

  // set up scenes
  for (auto sid = 0; sid < gltf->scenes_count; sid++) {
    auto  gscn = &gltf->scenes[sid];
    auto& scn  = scene.scenes.emplace_back();
    scn.name   = gscn->name ? gscn->name : "";
    for (auto nid = 0; nid < gscn->nodes_count; nid++) {
      scn.nodes.push_back(nmap.at(gscn->nodes[nid]));
    }
  }

  // update transforms
  for (auto& node : scene.nodes)
    if (node.parent < 0) update_transforms(scene, node);

#if 0
  // hasher for later
  struct sampler_map_hash {
    size_t operator()(
        const pair<cgltf_animation_sampler*, cgltf_animation_path_type>& value)
        const {
      auto hasher1 = std::hash<cgltf_animation_sampler*>();
      auto hasher2 = std::hash<int>();
      auto h       = (size_t)0;
      h ^= hasher1(value.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher2(value.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  // convert animations
  for (auto gid = 0; gid < gltf->animations_count; gid++) {
    auto ganm = &gltf->animations[gid];
    auto aid  = 0;
    auto sampler_map =
        unordered_map<pair<cgltf_animation_sampler*, cgltf_animation_path_type>,
            int, sampler_map_hash>();
    for (auto cid = 0; cid < ganm->channels_count; cid++) {
      auto gchannel = &ganm->channels[cid];
      auto path     = gchannel->target_path;
      if (sampler_map.find({gchannel->sampler, path}) == sampler_map.end()) {
        auto gsampler  = gchannel->sampler;
        auto animation = gltf_animation{};
        animation.uri  = (ganm->name ? ganm->name : "anim") +
                        std::to_string(aid++);
        animation.group = ganm->name ? ganm->name : "";
        auto input_view = accessor_values(gsampler->input);
        animation.times.resize(input_view.size());
        for (auto i = 0; i < input_view.size(); i++)
          animation.times[i] = input_view[i][0];
        switch (gsampler->interpolation) {
          case cgltf_interpolation_type_linear:
            animation.interpolation =
                gltf_animation::interpolation_type::linear;
            break;
          case cgltf_interpolation_type_step:
            animation.interpolation = gltf_animation::interpolation_type::step;
            break;
          case cgltf_interpolation_type_cubic_spline:
            animation.interpolation =
                gltf_animation::interpolation_type::bezier;
            break;
        }
        auto output_view = accessor_values(gsampler->output);
        switch (path) {
          case cgltf_animation_path_type_translation: {
            animation.translations.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.translations.push_back({(float)output_view[i][0],
                  (float)output_view[i][1], (float)output_view[i][2]});
          } break;
          case cgltf_animation_path_type_rotation: {
            animation.rotations.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.rotations.push_back(
                  {(float)output_view[i][0], (float)output_view[i][1],
                      (float)output_view[i][2], (float)output_view[i][3]});
          } break;
          case cgltf_animation_path_type_scale: {
            animation.scales.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.scales.push_back({(float)output_view[i][0],
                  (float)output_view[i][1], (float)output_view[i][2]});
          } break;
          case cgltf_animation_path_type_weights: {
            throw std::runtime_error("weights not supported for now");
                    // // get a node that it refers to
                    // auto ncomp = 0;
                    // auto gnode = gltf->get(gchannel->target->node);
                    // auto gmesh = gltf->get(gnode->mesh);
                    // if (gmesh) {
                    //     for (auto gshp : gmesh->primitives) {
                    //         ncomp = max((int)gshp->targets.size(), ncomp);
                    //     }
                    // }
                    // if (ncomp) {
                    //     auto values = vector<float>();
                    //     values.reserve(output_view.size());
                    //     for (auto i = 0; i < output_view.size(); i++)
                    //         values.push_back(output_view.get(i));
                    //     animation.weights.resize(values.size() / ncomp);
                    //     for (auto i = 0; i < animation.weights.size(); i++) {
                    //         animation.weights[i].resize(ncomp);
                    //         for (auto j = 0; j < ncomp; j++)
                    //             animation.weights[i][j] = values[i * ncomp + j];
                    //     }
                    // }
          } break;
          default: {
            throw std::runtime_error("bad gltf animation");
          }
        }
        sampler_map[{gchannel->sampler, path}] = (int)scene.animations.size();
        scene.animations.push_back(animation);
      }
      scene.animations[sampler_map.at({gchannel->sampler, path})]
          .targets.push_back(nmap.at(gchannel->target_node));
    }
  }
#endif
}

}  // namespace yocto

#if 0

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CYHAIR
// -----------------------------------------------------------------------------
namespace yocto {

struct cyhair_strand {
  vector<vec3f> positions;
  vector<float> radius;
  vector<float> transparency;
  vector<vec3f> color;
};

struct cyhair_data {
  vector<cyhair_strand> strands              = {};
  float                 default_thickness    = 0;
  float                 default_transparency = 0;
  vec3f                 default_color        = zero3f;
};

static void load_cyhair(const string& filename, cyhair_data& hair) {
  // open file
  hair     = {};
  auto fs_ = open_file(filename, "b");
  auto fs  = fs_.fs;

  // Bytes 0-3    Must be "HAIR" in ascii code (48 41 49 52)
  // Bytes 4-7    Number of hair strands as unsigned int
  // Bytes 8-11    Total number of points of all strands as unsigned int
  // Bytes 12-15    Bit array of data in the file
  // Bit-0 is 1 if the file has segments array.
  // Bit-1 is 1 if the file has points array (this bit must be 1).
  // Bit-2 is 1 if the file has radius array.
  // Bit-3 is 1 if the file has transparency array.
  // Bit-4 is 1 if the file has color array.
  // Bit-5 to Bit-31 are reserved for future extension (must be 0).
  // Bytes 16-19    Default number of segments of hair strands as unsigned int
  // If the file does not have a segments array, this default value is used.
  // Bytes 20-23    Default radius hair strands as float
  // If the file does not have a radius array, this default value is used.
  // Bytes 24-27    Default transparency hair strands as float
  // If the file does not have a transparency array, this default value is
  // used. Bytes 28-39    Default color hair strands as float array of size 3
  // If the file does not have a radius array, this default value is used.
  // Bytes 40-127    File information as char array of size 88 in ascii

  auto read_value = [](FILE* fs, auto& value) {
    if (fread(&value, sizeof(value), 1, fs) != 1) {
      throw std::runtime_error("cannot read from file");
    }
  };
  auto read_values = [](FILE* fs, auto& values) {
    if (values.empty()) return;
    if (fread(values.data(), sizeof(values[0]), values.size(), fs) !=
        values.size()) {
      throw std::runtime_error("cannot read from file");
    }
  };

  // parse header
  hair = cyhair_data{};
  struct cyhair_header {
    char         magic[4]             = {0};
    unsigned int num_strands          = 0;
    unsigned int num_points           = 0;
    unsigned int flags                = 0;
    unsigned int default_segments     = 0;
    float        default_thickness    = 0;
    float        default_transparency = 0;
    vec3f        default_color        = zero3f;
    char         info[88]             = {0};
  };
  static_assert(sizeof(cyhair_header) == 128);
  auto header = cyhair_header{};
  read_value(fs, header);
  if (header.magic[0] != 'H' || header.magic[1] != 'A' ||
      header.magic[2] != 'I' || header.magic[3] != 'R')
    throw std::runtime_error("bad cyhair header");

  // set up data
  hair.default_thickness    = header.default_thickness;
  hair.default_transparency = header.default_transparency;
  hair.default_color        = header.default_color;
  hair.strands.resize(header.num_strands);

  // get segments length
  auto segments = vector<unsigned short>();
  if (header.flags & 1) {
    segments.resize(header.num_strands);
    read_values(fs, segments);
  } else {
    segments.assign(header.num_strands, header.default_segments);
  }

  // check segment length
  auto total_length = 0;
  for (auto segment : segments) total_length += segment + 1;
  if (total_length != header.num_points) {
    throw std::runtime_error("bad cyhair file");
  }

  // read positions data
  if (header.flags & 2) {
    for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
      auto strand_size = (int)segments[strand_id] + 1;
      hair.strands[strand_id].positions.resize(strand_size);
      read_values(fs, hair.strands[strand_id].positions);
    }
  }
  // read radius data
  if (header.flags & 4) {
    for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
      auto strand_size = (int)segments[strand_id] + 1;
      hair.strands[strand_id].radius.resize(strand_size);
      read_values(fs, hair.strands[strand_id].radius);
    }
  }
  // read transparency data
  if (header.flags & 8) {
    for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
      auto strand_size = (int)segments[strand_id] + 1;
      hair.strands[strand_id].transparency.resize(strand_size);
      read_values(fs, hair.strands[strand_id].transparency);
    }
  }
  // read color data
  if (header.flags & 16) {
    for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
      auto strand_size = (int)segments[strand_id] + 1;
      hair.strands[strand_id].color.resize(strand_size);
      read_values(fs, hair.strands[strand_id].color);
    }
  }
}

// Compute per-vertex tangents for lines.
static void compute_cyhair_tangents(vector<vec3f>& tangents,
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  if (tangents.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& tangent : tangents) tangent = zero3f;
  for (auto& l : lines) {
    auto tangent = normalize(positions[l.y] - positions[l.x]);
    auto len     = length(positions[l.y] - positions[l.x]);
    tangents[l.x] += tangent * len;
    tangents[l.y] += tangent * len;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
}

void load_cyhair_shape(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec4f>& color, vector<float>& radius, bool flip_texcoord) {
  // load hair file
  auto hair = cyhair_data();
  load_cyhair(filename, hair);

  // generate curve data
  for (auto& strand : hair.strands) {
    auto offset = (int)positions.size();
    for (auto segment = 0; segment < (int)strand.positions.size() - 1;
         segment++) {
      lines.push_back({offset + segment, offset + segment + 1});
    }
    positions.insert(
        positions.end(), strand.positions.begin(), strand.positions.end());
    if (strand.radius.empty()) {
      radius.insert(
          radius.end(), strand.positions.size(), hair.default_thickness);
    } else {
      radius.insert(radius.end(), strand.radius.begin(), strand.radius.end());
    }
    if (strand.color.empty()) {
      color.insert(color.end(), strand.positions.size(),
          {hair.default_color.x, hair.default_color.y, hair.default_color.z,
              1});
    } else {
      for (auto i = 0; i < strand.color.size(); i++) {
        auto scolor = strand.color[i];
        color.push_back({scolor.x, scolor.y, scolor.z, 1});
      }
    }
  }

  // flip yz
  for (auto& p : positions) std::swap(p.y, p.z);

  // compute tangents
  normals.resize(positions.size());
  compute_cyhair_tangents(normals, lines, positions);

  // fix colors
  for (auto& c : color) c = {pow(xyz(c), 2.2f), c.w};
}

}  // namespace yocto

#endif
