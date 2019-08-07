//
// Implementation for Yocto/Pbrt
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

#include "yocto_pbrt.h"
#include "yocto_image.h"

#include <ctype.h>
#include <string_view>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF LOW LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Type aliases for readability
using string_view = std::string_view;
using namespace std::literals::string_view_literals;

// A file holder that closes a file when destructed. Useful for RIIA
struct file_holder {
  FILE*  fs       = nullptr;
  string filename = "";

  file_holder() {}
  file_holder(file_holder&& other) {
    this->fs       = other.fs;
    this->filename = other.filename;
    other.fs       = nullptr;
  }
  file_holder(const file_holder&) = delete;
  file_holder& operator=(const file_holder&) = delete;
  ~file_holder() {
    if (fs) fclose(fs);
  }
};

static inline void open_input_file(
    file_holder& file, const string& filename, bool binary = false) {
  auto& fs = file.fs;
  fs       = fopen(filename.c_str(), !binary ? "rt" : "rb");
  if (!fs) throw std::runtime_error("could not open file " + filename);
}

// Read a line
static inline bool read_line(FILE* fs, char* buffer, size_t size) {
  return fgets(buffer, size, fs) != nullptr;
}

static inline bool is_pbrt_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_pbrt_newline(char c) { return c == '\r' || c == '\n'; }

static inline void skip_pbrt_whitespace(string_view& str) {
  while (!str.empty() && is_pbrt_space(str.front())) str.remove_prefix(1);
}
static inline void remove_pbrt_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_pbrt_newline(str.back())) str.remove_suffix(1);
  auto cpy       = str;
  auto in_string = false;
  while (!cpy.empty()) {
    if (cpy.front() == '"') in_string = !in_string;
    if (cpy.front() == comment_char && !in_string) break;
    cpy.remove_prefix(1);
  }
  str.remove_suffix(cpy.size());
}

// Read a pbrt command from file
bool read_pbrt_line(FILE* fs, string& cmd) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs);
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_pbrt_comment(line);
    skip_pbrt_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs, pos, SEEK_SET);
        return true;
      } else {
        found = true;
      }
    } else if (!found) {
      throw std::runtime_error("bad pbrt command");
    }
    cmd += line;
    cmd += " ";
    pos = ftell(fs);
  }
  return found;
}

// parse a quoted string
static inline void parse_pbrt_value(string_view& str, string_view& value) {
  skip_pbrt_whitespace(str);
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

static inline void parse_pbrt_value(string_view& str, string& value) {
  auto view = ""sv;
  parse_pbrt_value(str, view);
  value = string{view};
}

// parse a quoted string
static inline void parse_pbrt_command(string_view& str, string& value) {
  skip_pbrt_whitespace(str);
  if (!isalpha((int)str.front())) {
    throw std::runtime_error("bad command");
  }
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == string_view::npos) {
    value.assign(str);
    str.remove_prefix(str.size());
  } else {
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
  }
}

// parse a number
static inline void parse_pbrt_value(string_view& str, float& value) {
  skip_pbrt_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtof(str.data(), &next);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}

// parse a number
static inline void parse_pbrt_value(string_view& str, int& value) {
  skip_pbrt_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtol(str.data(), &next, 10);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}
static inline void parse_pbrt_value(string_view& str, bool& value) {
  auto value_name = ""s;
  parse_pbrt_value(str, value_name);
  if (value_name == "true") {
    value = true;
  } else if (value_name == "false") {
    value = false;
  } else {
    throw std::runtime_error("expected boolean");
  }
}
template <typename T>
static inline void parse_pbrt_value(
    string_view& str, T& value, unordered_map<string, T>& value_names) {
  auto value_name = ""s;
  parse_pbrt_value(str, value_name);
  try {
    value = value_names.at(value_name);
  } catch (std::out_of_range&) {
    throw std::runtime_error("expected enum value");
  }
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::bilerp_t::mapping_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::bilerp_t::mapping_type>{
          {"uv", pbrt_texture::bilerp_t::mapping_type::uv},
          {"spherical", pbrt_texture::bilerp_t::mapping_type::spherical},
          {"cylindrical", pbrt_texture::bilerp_t::mapping_type::cylindrical},
          {"planar", pbrt_texture::bilerp_t::mapping_type::planar},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::checkerboard_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::dots_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::imagemap_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::uv_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}

static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::checkerboard_t::aamode_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::checkerboard_t::aamode_type>{
          {"closedform", pbrt_texture::checkerboard_t::aamode_type::closedform},
          {"none", pbrt_texture::checkerboard_t::aamode_type::none},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::imagemap_t::wrap_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::imagemap_t::wrap_type>{
          {"repeat", pbrt_texture::imagemap_t::wrap_type::repeat},
          {"clamp", pbrt_texture::imagemap_t::wrap_type::clamp},
          {"black", pbrt_texture::imagemap_t::wrap_type::black},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_shape::curve_t::basis_t& value) {
  static auto value_names = unordered_map<string, pbrt_shape::curve_t::basis_t>{
      {"bezier", pbrt_shape::curve_t::basis_t::bezier},
      {"bspline", pbrt_shape::curve_t::basis_t::bspline},
  };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_shape::curve_t::type_t& value) {
  static auto value_names = unordered_map<string, pbrt_shape::curve_t::type_t>{
      {"flat", pbrt_shape::curve_t::type_t::flat},
      {"cylinder", pbrt_shape::curve_t::type_t::cylinder},
      {"ribbon", pbrt_shape::curve_t::type_t::ribbon},
  };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_accelerator::bvh_t::splitmethod_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_accelerator::bvh_t::splitmethod_t>{
          {"sah", pbrt_accelerator::bvh_t::splitmethod_t::sah},
          {"equal", pbrt_accelerator::bvh_t::splitmethod_t::equal},
          {"middle", pbrt_accelerator::bvh_t::splitmethod_t::middle},
          {"hlbvh", pbrt_accelerator::bvh_t::splitmethod_t::hlbvh},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::path_t::lightsamplestrategy_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_integrator::path_t::lightsamplestrategy_t>{
          {"power", pbrt_integrator::path_t::lightsamplestrategy_t::power},
          {"spatial", pbrt_integrator::path_t::lightsamplestrategy_t::spatial},
          {"uniform", pbrt_integrator::path_t::lightsamplestrategy_t::uniform},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(string_view&       str,
    pbrt_integrator::volpath_t::lightsamplestrategy_t& value) {
  return parse_pbrt_value(
      str, (pbrt_integrator::path_t::lightsamplestrategy_t&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::bdpt_t::lightsamplestrategy_t& value) {
  return parse_pbrt_value(
      str, (pbrt_integrator::path_t::lightsamplestrategy_t&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::directlighting_t::strategy_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_integrator::directlighting_t::strategy_t>{
          {"all", pbrt_integrator::directlighting_t::strategy_t::all},
          {"one", pbrt_integrator::directlighting_t::strategy_t::one},
      };
  return parse_pbrt_value(str, value, value_names);
}

// parse a vec type
static inline void parse_pbrt_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec3i& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec4i& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, pbrt_spectrum3f& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}

// Check next
static inline bool is_pbrt_string(string_view& str) {
  skip_pbrt_whitespace(str);
  return !str.empty() && str.front() == '"';
}
static inline bool is_open_bracket(string_view& str) {
  skip_pbrt_whitespace(str);
  return !str.empty() && str.front() == '[';
}
static inline bool is_close_bracket(string_view& str) {
  skip_pbrt_whitespace(str);
  return !str.empty() && str.front() == ']';
}
static inline bool is_pbrt_param(string_view& str) {
  skip_pbrt_whitespace(str);
  return is_pbrt_string(str);
}

// parse a quoted string
static inline void parse_pbrt_nametype(
    string_view& str_, string& name, string& type) {
  auto value = ""s;
  parse_pbrt_value(str_, value);
  auto str  = string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == string_view::npos) {
    throw std::runtime_error("bad type " + value);
  }
  type = string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == string_view::npos) {
    throw std::runtime_error("bad type " + value);
  }
  str.remove_prefix(pos2);
  name = string(str);
}

static inline void skip_pbrt_open_bracket(string_view& str) {
  if (!is_open_bracket(str)) throw std::runtime_error("expected bracket");
  str.remove_prefix(1);
  skip_pbrt_whitespace(str);
}
static inline void skip_pbrt_close_bracket(string_view& str) {
  if (!is_close_bracket(str)) throw std::runtime_error("expected bracket");
  str.remove_prefix(1);
  skip_pbrt_whitespace(str);
}

template <typename T>
static inline void parse_pbrt_param(string_view& str, T& value) {
  auto has_brackets = is_open_bracket(str);
  if (has_brackets) skip_pbrt_open_bracket(str);
  parse_pbrt_value(str, value);
  if (has_brackets) skip_pbrt_close_bracket(str);
}

template <typename T>
static inline void parse_pbrt_param(string_view& str, vector<T>& values) {
  skip_pbrt_open_bracket(str);
  values.clear();
  while (!is_close_bracket(str)) {
    values.push_back({});
    parse_pbrt_value(str, values.back());
  }
  skip_pbrt_close_bracket(str);
}

template <typename T>
static inline bool is_pbrt_type_compatible(const string& type) {
  if constexpr (std::is_same<T, int>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, float>::value) {
    return type == "float";
  } else if constexpr (std::is_same<T, bool>::value) {
    return type == "bool";
  } else if constexpr (std::is_same<T, string>::value) {
    return type == "string";
  } else if constexpr (std::is_same<T, vec2f>::value) {
    return type == "point2" || type == "vector2" || type == "float";
  } else if constexpr (std::is_same<T, vec3f>::value) {
    return type == "point3" || type == "vector3" || type == "normal3" ||
           type == "point" || type == "vector" || type == "normal" ||
           type == "float";
  } else if constexpr (std::is_same<T, vec4f>::value) {
    return type == "float";
  } else if constexpr (std::is_same<T, pbrt_spectrum3f>::value) {
    return type == "rgb" || type == "pbrt_spectrum" || type == "blackbody";
  } else if constexpr (std::is_same<T, vec3i>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, vec4i>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, bbox2f>::value) {
    return type == "float";
  } else if constexpr (std::is_enum<T>::value) {
    return type == "string";
  } else {
    return false;
  }
}

template <typename T>
static inline void parse_pbrt_param(
    string_view& str, const string& type, T& value) {
  if (!is_pbrt_type_compatible<T>(type)) {
    throw std::runtime_error("incompatible type " + type);
  }
  parse_pbrt_param(str, value);
}

static inline pair<vec3f, vec3f> get_pbrt_etak(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> metal_ior_table = {
      {"a-C", {{2.9440999183f, 2.2271502925f, 1.9681668794f},
                  {0.8874329109f, 0.7993216383f, 0.8152862927f}}},
      {"Ag", {{0.1552646489f, 0.1167232965f, 0.1383806959f},
                 {4.8283433224f, 3.1222459278f, 2.1469504455f}}},
      {"Al", {{1.6574599595f, 0.8803689579f, 0.5212287346f},
                 {9.2238691996f, 6.2695232477f, 4.8370012281f}}},
      {"AlAs", {{3.6051023902f, 3.2329365777f, 2.2175611545f},
                   {0.0006670247f, -0.0004999400f, 0.0074261204f}}},
      {"AlSb", {{-0.0485225705f, 4.1427547893f, 4.6697691348f},
                   {-0.0363741915f, 0.0937665154f, 1.3007390124f}}},
      {"Au", {{0.1431189557f, 0.3749570432f, 1.4424785571f},
                 {3.9831604247f, 2.3857207478f, 1.6032152899f}}},
      {"Be", {{4.1850592788f, 3.1850604423f, 2.7840913457f},
                 {3.8354398268f, 3.0101260162f, 2.8690088743f}}},
      {"Cr", {{4.3696828663f, 2.9167024892f, 1.6547005413f},
                 {5.2064337956f, 4.2313645277f, 3.7549467933f}}},
      {"CsI", {{2.1449030413f, 1.7023164587f, 1.6624194173f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Cu", {{0.2004376970f, 0.9240334304f, 1.1022119527f},
                 {3.9129485033f, 2.4528477015f, 2.1421879552f}}},
      {"Cu2O", {{3.5492833755f, 2.9520622449f, 2.7369202137f},
                   {0.1132179294f, 0.1946659670f, 0.6001681264f}}},
      {"CuO", {{3.2453822204f, 2.4496293965f, 2.1974114493f},
                  {0.5202739621f, 0.5707372756f, 0.7172250613f}}},
      {"d-C", {{2.7112524747f, 2.3185812849f, 2.2288565009f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Hg", {{2.3989314904f, 1.4400254917f, 0.9095512090f},
                 {6.3276269444f, 4.3719414152f, 3.4217899270f}}},
      {"HgTe", {{4.7795267752f, 3.2309984581f, 2.6600252401f},
                   {1.6319827058f, 1.5808189339f, 1.7295753852f}}},
      {"Ir", {{3.0864098394f, 2.0821938440f, 1.6178866805f},
                 {5.5921510077f, 4.0671757150f, 3.2672611269f}}},
      {"K", {{0.0640493070f, 0.0464100621f, 0.0381842017f},
                {2.1042155920f, 1.3489364357f, 0.9132113889f}}},
      {"Li", {{0.2657871942f, 0.1956102432f, 0.2209198538f},
                 {3.5401743407f, 2.3111306542f, 1.6685930000f}}},
      {"MgO", {{2.0895885542f, 1.6507224525f, 1.5948759692f},
                  {0.0000000000f, -0.0000000000f, 0.0000000000f}}},
      {"Mo", {{4.4837010280f, 3.5254578255f, 2.7760769438f},
                 {4.1111307988f, 3.4208716252f, 3.1506031404f}}},
      {"Na", {{0.0602665320f, 0.0561412435f, 0.0619909494f},
                 {3.1792906496f, 2.1124800781f, 1.5790940266f}}},
      {"Nb", {{3.4201353595f, 2.7901921379f, 2.3955856658f},
                 {3.4413817900f, 2.7376437930f, 2.5799132708f}}},
      {"Ni", {{2.3672753521f, 1.6633583302f, 1.4670554172f},
                 {4.4988329911f, 3.0501643957f, 2.3454274399f}}},
      {"Rh", {{2.5857954933f, 1.8601866068f, 1.5544279524f},
                 {6.7822927110f, 4.7029501026f, 3.9760892461f}}},
      {"Se-e", {{5.7242724833f, 4.1653992967f, 4.0816099264f},
                   {0.8713747439f, 1.1052845009f, 1.5647788766f}}},
      {"Se", {{4.0592611085f, 2.8426947380f, 2.8207582835f},
                 {0.7543791750f, 0.6385150558f, 0.5215872029f}}},
      {"SiC", {{3.1723450205f, 2.5259677964f, 2.4793623897f},
                  {0.0000007284f, -0.0000006859f, 0.0000100150f}}},
      {"SnTe", {{4.5251865890f, 1.9811525984f, 1.2816819226f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Ta", {{2.0625846607f, 2.3930915569f, 2.6280684948f},
                 {2.4080467973f, 1.7413705864f, 1.9470377016f}}},
      {"Te-e", {{7.5090397678f, 4.2964603080f, 2.3698732430f},
                   {5.5842076830f, 4.9476231084f, 3.9975145063f}}},
      {"Te", {{7.3908396088f, 4.4821028985f, 2.6370708478f},
                 {3.2561412892f, 3.5273908133f, 3.2921683116f}}},
      {"ThF4", {{1.8307187117f, 1.4422274283f, 1.3876488528f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"TiC", {{3.7004673762f, 2.8374356509f, 2.5823030278f},
                  {3.2656905818f, 2.3515586388f, 2.1727857800f}}},
      {"TiN", {{1.6484691607f, 1.1504482522f, 1.3797795097f},
                  {3.3684596226f, 1.9434888540f, 1.1020123347f}}},
      {"TiO2-e", {{3.1065574823f, 2.5131551146f, 2.5823844157f},
                     {0.0000289537f, -0.0000251484f, 0.0001775555f}}},
      {"TiO2", {{3.4566203131f, 2.8017076558f, 2.9051485020f},
                   {0.0001026662f, -0.0000897534f, 0.0006356902f}}},
      {"VC", {{3.6575665991f, 2.7527298065f, 2.5326814570f},
                 {3.0683516659f, 2.1986687713f, 1.9631816252f}}},
      {"VN", {{2.8656011588f, 2.1191817791f, 1.9400767149f},
                 {3.0323264950f, 2.0561075580f, 1.6162930914f}}},
      {"V", {{4.2775126218f, 3.5131538236f, 2.7611257461f},
                {3.4911844504f, 2.8893580874f, 3.1116965117f}}},
      {"W", {{4.3707029924f, 3.3002972445f, 2.9982666528f},
                {3.5006778591f, 2.6048652781f, 2.2731930614f}}},
  };
  return metal_ior_table.at(name);
}

static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_spectrum3f& value) {
  bool verbose = false;
  if (type == "rgb") {
    parse_pbrt_param(str, value);
  } else if (type == "color") {
    parse_pbrt_param(str, value);
  } else if (type == "float") {
    auto valuef = 0.0f;
    parse_pbrt_param(str, valuef);
    value = {valuef, valuef, valuef};
  } else if (type == "blackbody") {
    auto blackbody = zero2f;
    parse_pbrt_param(str, blackbody);
    (vec3f&)value = blackbody_to_rgb(blackbody.x) * blackbody.y;
  } else if (type == "spectrum" && is_pbrt_string(str)) {
    if (verbose) printf("spectrum  not well supported\n");
    auto filename = ""s;
    parse_pbrt_param(str, filename);
    auto filenamep = fs::path(filename).filename();
    if (filenamep.extension() == ".spd") {
      filenamep = filenamep.replace_extension("");
      if (filenamep == "SHPS") {
        value = {1, 1, 1};
      } else if (filenamep.extension() == ".eta") {
        auto eta = get_pbrt_etak(filenamep.replace_extension("")).first;
        value    = {eta.x, eta.y, eta.z};
      } else if (filenamep.extension() == ".k") {
        auto k = get_pbrt_etak(filenamep.replace_extension("")).second;
        value  = {k.x, k.y, k.z};
      } else {
        throw std::runtime_error("unknown spectrum file " + filename);
      }
    } else {
      throw std::runtime_error("unsupported spectrum format");
      // value = {1, 0, 0};
    }
  } else if (type == "spectrum" && !is_pbrt_string(str)) {
    if (verbose) printf("spectrum  not well supported\n");
    auto values = vector<float>{};
    parse_pbrt_param(str, values);
    value = {1, 0, 0};
  } else {
    throw std::runtime_error("unsupported spectrum type");
  }
}

template <typename T>
static inline void parse_pbrt_param(
    string_view& str, const string& type, vector<T>& value) {
  if (!is_pbrt_type_compatible<T>(type)) {
    throw std::runtime_error("incompatible type " + type);
  }
  parse_pbrt_param(str, value);
}

static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_textured3f& value) {
  if (type == "texture") {
    parse_pbrt_param(str, value.texture);
  } else {
    parse_pbrt_param(str, type, value.value);
  }
}
static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_textured1f& value) {
  if (type == "texture") {
    parse_pbrt_param(str, value.texture);
  } else {
    parse_pbrt_param(str, type, value.value);
  }
}

static inline void skip_pbrt_value(string_view& str) {
  skip_pbrt_whitespace(str);
  if (str.front() == '"') {
    str.remove_prefix(1);
    str.remove_prefix(str.find('"') + 1);
  } else {
    str.remove_prefix(str.find_first_of(" \n\t\r],\""));
  }
  skip_pbrt_whitespace(str);
}

static inline void skip_pbrt_param(string_view& str) {
  if (is_open_bracket(str)) {
    skip_pbrt_open_bracket(str);
    while (!is_close_bracket(str)) skip_pbrt_value(str);
    skip_pbrt_close_bracket(str);
  } else {
    skip_pbrt_value(str);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Parse Accelerator
static inline void parse_pbrt_accelerator(
    string_view& str, const string& type, pbrt_accelerator& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "bvh") {
    auto tvalue = pbrt_accelerator::bvh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxnodeprims") {
        parse_pbrt_param(str, ptype, tvalue.maxnodeprims);
      } else if (pname == "splitmethod") {
        parse_pbrt_param(str, ptype, tvalue.splitmethod);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_accelerator::type_t::bvh;
    value.bvh  = tvalue;
  } else if (type == "kdtree") {
    auto tvalue = pbrt_accelerator::kdtree_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "intersectcost") {
        parse_pbrt_param(str, ptype, tvalue.intersectcost);
      } else if (pname == "traversalcost") {
        parse_pbrt_param(str, ptype, tvalue.traversalcost);
      } else if (pname == "emptybonus") {
        parse_pbrt_param(str, ptype, tvalue.emptybonus);
      } else if (pname == "maxprims") {
        parse_pbrt_param(str, ptype, tvalue.maxprims);
      } else if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_accelerator::type_t::kdtree;
    value.kdtree = tvalue;
  } else {
    throw std::runtime_error("unknown Accelerator " + type);
  }
}

// Parse Integrator
static inline void parse_pbrt_integrator(
    string_view& str, const string& type, pbrt_integrator& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "path") {
    auto tvalue = pbrt_integrator::path_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "rrthreshold") {
        parse_pbrt_param(str, ptype, tvalue.rrthreshold);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
      // parse_pbrt_optional_param(str, "lightsamplestrategy",
      // tvalue.lightsamplestrategy); // TODO: enums
    }
    value.type = pbrt_integrator::type_t::path;
    value.path = tvalue;
  } else if (type == "volpath") {
    auto tvalue = pbrt_integrator::volpath_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "rrthreshold") {
        parse_pbrt_param(str, ptype, tvalue.rrthreshold);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_integrator::type_t::volpath;
    value.volpath = tvalue;
  } else if (type == "directlighting") {
    auto tvalue = pbrt_integrator::directlighting_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "strategy") {
        parse_pbrt_param(str, ptype, tvalue.strategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type           = pbrt_integrator::type_t::directlighting;
    value.directlighting = tvalue;
  } else if (type == "bdpt") {
    auto tvalue = pbrt_integrator::bdpt_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else if (pname == "visualizestrategies") {
        parse_pbrt_param(str, ptype, tvalue.visualizestrategies);
      } else if (pname == "visualizeweights") {
        parse_pbrt_param(str, ptype, tvalue.visualizeweights);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::bdpt;
    value.bdpt = tvalue;
  } else if (type == "mlt") {
    auto tvalue = pbrt_integrator::mlt_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "bootstrapsamples") {
        parse_pbrt_param(str, ptype, tvalue.bootstrapsamples);
      } else if (pname == "chains") {
        parse_pbrt_param(str, ptype, tvalue.chains);
      } else if (pname == "mutationsperpixel") {
        parse_pbrt_param(str, ptype, tvalue.mutationsperpixel);
      } else if (pname == "largestepprobability") {
        parse_pbrt_param(str, ptype, tvalue.largestepprobability);
      } else if (pname == "sigma") {
        parse_pbrt_param(str, ptype, tvalue.sigma);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::mlt;
    value.mlt  = tvalue;
  } else if (type == "sppm") {
    auto tvalue = pbrt_integrator::sppm_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "iterations") {
        parse_pbrt_param(str, ptype, tvalue.iterations);
      } else if (pname == "numiterations") {
        parse_pbrt_param(str, ptype, tvalue.iterations);
      } else if (pname == "photonsperiteration") {
        parse_pbrt_param(str, ptype, tvalue.photonsperiteration);
      } else if (pname == "imagewritefrequency") {
        parse_pbrt_param(str, ptype, tvalue.imagewritefrequency);
      } else if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::sppm;
    value.sppm = tvalue;
  } else if (type == "whitted") {
    auto tvalue = pbrt_integrator::whitted_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_integrator::type_t::whitted;
    value.whitted = tvalue;
  } else {
    throw std::runtime_error("unknown Integrator " + type);
  }
}

// Parse Sampler
static inline void parse_pbrt_sampler(
    string_view& str, const string& type, pbrt_sampler& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "random") {
    auto tvalue = pbrt_sampler::random_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_sampler::type_t::random;
    value.random = tvalue;
  } else if (type == "halton") {
    auto tvalue = pbrt_sampler::halton_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_sampler::type_t::halton;
    value.halton = tvalue;
  } else if (type == "sobol") {
    auto tvalue = pbrt_sampler::sobol_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_sampler::type_t::sobol;
    value.sobol = tvalue;
  } else if (type == "02sequence") {
    auto tvalue = pbrt_sampler::zerotwosequence_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type            = pbrt_sampler::type_t::zerotwosequence;
    value.zerotwosequence = tvalue;
  } else if (type == "lowdiscrepancy") {
    auto tvalue = pbrt_sampler::zerotwosequence_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type            = pbrt_sampler::type_t::zerotwosequence;
    value.zerotwosequence = tvalue;
  } else if (type == "maxmindist") {
    auto tvalue = pbrt_sampler::maxmindist_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_sampler::type_t::maxmindist;
    value.maxmindist = tvalue;
  } else if (type == "stratified") {
    auto tvalue = pbrt_sampler::stratified_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xsamples") {
        parse_pbrt_param(str, ptype, tvalue.xsamples);
      } else if (pname == "ysamples") {
        parse_pbrt_param(str, ptype, tvalue.ysamples);
      } else if (pname == "jitter") {
        parse_pbrt_param(str, ptype, tvalue.jitter);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_sampler::type_t::stratified;
    value.stratified = tvalue;
  } else {
    throw std::runtime_error("unknown Sampler " + type);
  }
}

// Parse Filter
static inline void parse_pbrt_filter(
    string_view& str, const string& type, pbrt_filter& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "box") {
    auto tvalue = pbrt_filter::box_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_filter::type_t::box;
    value.box  = tvalue;
  } else if (type == "gaussian") {
    auto tvalue = pbrt_filter::gaussian_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::gaussian;
    value.gaussian = tvalue;
  } else if (type == "mitchell") {
    auto tvalue = pbrt_filter::mitchell_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "B") {
        parse_pbrt_param(str, ptype, tvalue.B);
      } else if (pname == "C") {
        parse_pbrt_param(str, ptype, tvalue.C);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::mitchell;
    value.mitchell = tvalue;
  } else if (type == "sinc") {
    auto tvalue = pbrt_filter::sinc_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "tau") {
        parse_pbrt_param(str, ptype, tvalue.tau);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_filter::type_t::sinc;
    value.sinc = tvalue;
  } else if (type == "triangle") {
    auto tvalue = pbrt_filter::triangle_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::triangle;
    value.triangle = tvalue;
  } else {
    throw std::runtime_error("unknown PixelFilter " + type);
  }
}

// Parse Filter
static inline void parse_pbrt_film(
    string_view& str, const string& type, pbrt_film& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "image") {
    auto tvalue = pbrt_film::image_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xresolution") {
        parse_pbrt_param(str, ptype, tvalue.xresolution);
      } else if (pname == "yresolution") {
        parse_pbrt_param(str, ptype, tvalue.yresolution);
      } else if (pname == "yresolution") {
        parse_pbrt_param(str, ptype, tvalue.yresolution);
      } else if (pname == "cropwindow") {
        parse_pbrt_param(str, ptype, tvalue.cropwindow);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "maxsampleluminance") {
        parse_pbrt_param(str, ptype, tvalue.maxsampleluminance);
      } else if (pname == "diagonal") {
        parse_pbrt_param(str, ptype, tvalue.diagonal);
      } else if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_film::type_t::image;
    value.image = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse Camera
static inline void parse_pbrt_camera(
    string_view& str, const string& type, pbrt_camera& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "perspective") {
    auto tvalue = pbrt_camera::perspective_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "fov") {
        parse_pbrt_param(str, ptype, tvalue.fov);
      } else if (pname == "frameaspectratio") {
        parse_pbrt_param(str, ptype, tvalue.frameaspectratio);
      } else if (pname == "lensradius") {
        parse_pbrt_param(str, ptype, tvalue.lensradius);
      } else if (pname == "focaldistance") {
        parse_pbrt_param(str, ptype, tvalue.focaldistance);
      } else if (pname == "screenwindow") {
        parse_pbrt_param(str, ptype, tvalue.screenwindow);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_camera::type_t::perspective;
    value.perspective = tvalue;
  } else if (type == "orthographic") {
    auto tvalue = pbrt_camera::orthographic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "frameaspectratio") {
        parse_pbrt_param(str, ptype, tvalue.frameaspectratio);
      } else if (pname == "lensradius") {
        parse_pbrt_param(str, ptype, tvalue.lensradius);
      } else if (pname == "focaldistance") {
        parse_pbrt_param(str, ptype, tvalue.focaldistance);
      } else if (pname == "screenwindow") {
        parse_pbrt_param(str, ptype, tvalue.screenwindow);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_camera::type_t::orthographic;
    value.orthographic = tvalue;
  } else if (type == "environment") {
    auto tvalue = pbrt_camera::environment_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_camera::type_t::environment;
    value.environment = tvalue;
  } else if (type == "realistic") {
    auto tvalue = pbrt_camera::realistic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "lensfile") {
        parse_pbrt_param(str, ptype, tvalue.lensfile);
        // example: wide.22mm.dat
        auto lensfile = fs::path(tvalue.lensfile).filename().string();
        lensfile      = lensfile.substr(0, lensfile.size() - 4);
        lensfile      = lensfile.substr(lensfile.find('.') + 1);
        lensfile      = lensfile.substr(0, lensfile.size() - 2);
        tvalue.approx_focallength = std::atof(lensfile.c_str());
      } else if (pname == "aperturediameter") {
        parse_pbrt_param(str, ptype, tvalue.aperturediameter);
      } else if (pname == "focusdistance") {
        parse_pbrt_param(str, ptype, tvalue.focusdistance);
      } else if (pname == "simpleweighting") {
        parse_pbrt_param(str, ptype, tvalue.simpleweighting);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type      = pbrt_camera::type_t::realistic;
    value.realistic = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse Texture
static inline void parse_pbrt_texture(
    string_view& str, const string& type, pbrt_texture& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "constant") {
    auto tvalue = pbrt_texture::constant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "value") {
        parse_pbrt_param(str, ptype, tvalue.value);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::constant;
    value.constant = tvalue;
  } else if (type == "bilerp") {
    auto tvalue = pbrt_texture::bilerp_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "v00") {
        parse_pbrt_param(str, ptype, tvalue.v00);
      } else if (pname == "v01") {
        parse_pbrt_param(str, ptype, tvalue.v01);
      } else if (pname == "v10") {
        parse_pbrt_param(str, ptype, tvalue.v10);
      } else if (pname == "v11") {
        parse_pbrt_param(str, ptype, tvalue.v11);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_texture::type_t::bilerp;
    value.bilerp = tvalue;
  } else if (type == "checkerboard") {
    auto tvalue = pbrt_texture::checkerboard_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "dimension") {
        parse_pbrt_param(str, ptype, tvalue.dimension);
      } else if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else if (pname == "aamode") {
        parse_pbrt_param(str, ptype, tvalue.aamode);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_texture::type_t::checkerboard;
    value.checkerboard = tvalue;
  } else if (type == "dots") {
    auto tvalue = pbrt_texture::dots_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "inside") {
        parse_pbrt_param(str, ptype, tvalue.inside);
      } else if (pname == "outside") {
        parse_pbrt_param(str, ptype, tvalue.outside);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::dots;
    value.dots = tvalue;
  } else if (type == "imagemap") {
    auto tvalue = pbrt_texture::imagemap_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else if (pname == "wrap") {
        parse_pbrt_param(str, ptype, tvalue.wrap);
      } else if (pname == "maxanisotropy") {
        parse_pbrt_param(str, ptype, tvalue.maxanisotropy);
      } else if (pname == "trilinear") {
        parse_pbrt_param(str, ptype, tvalue.trilinear);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "gamma") {
        parse_pbrt_param(str, ptype, tvalue.gamma);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::imagemap;
    value.imagemap = tvalue;
  } else if (type == "mix") {
    auto tvalue = pbrt_texture::mix_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else if (pname == "amount") {
        parse_pbrt_param(str, ptype, tvalue.amount);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::mix;
    value.mix  = tvalue;
  } else if (type == "scale") {
    auto tvalue = pbrt_texture::scale_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_texture::type_t::scale;
    value.scale = tvalue;
  } else if (type == "fbm") {
    auto tvalue = pbrt_texture::fbm_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::fbm;
    value.fbm  = tvalue;
  } else if (type == "wrinkled") {
    auto tvalue = pbrt_texture::wrinkled_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::wrinkled;
    value.wrinkled = tvalue;
  } else if (type == "windy") {
    auto tvalue = pbrt_texture::windy_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "") {
        // TODO: missing params
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_texture::type_t::windy;
    value.windy = tvalue;
  } else if (type == "marble") {
    auto tvalue = pbrt_texture::marble_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "variation") {
        parse_pbrt_param(str, ptype, tvalue.variation);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_texture::type_t::marble;
    value.marble = tvalue;
  } else if (type == "uv") {
    auto tvalue = pbrt_texture::uv_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::uv;
    value.uv   = tvalue;
  } else {
    throw std::runtime_error("unknown Texture " + type);
  }
}

void approximate_fourier_material(pbrt_material::fourier_t& fourier) {
  auto filename = fs::path(fourier.bsdffile).filename().string();
  if (filename == "paint.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.6f, 0.6f};
    // plastic.Ks = {0.4f, 0.4f, 0.4f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.2f;
    plastic.vroughness = 0.2f;
  } else if (filename == "ceramic.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.6f, 0.6f};
    // plastic.Ks = {0.1f, 0.1f, 0.1f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.025f;
    plastic.vroughness = 0.025f;
  } else if (filename == "leather.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.57f, 0.48f};
    // plastic.Ks = {0.1f, 0.1f, 0.1f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.3f;
    plastic.vroughness = 0.3f;
  } else if (filename == "coated_copper.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::metal;
    auto& metal         = fourier.approx_metal;
    auto  etak          = get_pbrt_etak("Cu");
    metal.eta           = {etak.first.x, etak.first.y, etak.first.z};
    metal.k             = {etak.second.x, etak.second.y, etak.second.z};
    metal.uroughness    = 0.01f;
    metal.vroughness    = 0.01f;
  } else if (filename == "roughglass_alpha_0.2.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::glass;
    auto& glass         = fourier.approx_glass;
    glass.uroughness    = 0.2f;
    glass.vroughness    = 0.2f;
    glass.Kr            = {1, 1, 1};
    glass.Kt            = {1, 1, 1};
  } else if (filename == "roughgold_alpha_0.2.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::metal;
    auto& metal         = fourier.approx_metal;
    auto  etak          = get_pbrt_etak("Au");
    metal.eta           = {etak.first.x, etak.first.y, etak.first.z};
    metal.k             = {etak.second.x, etak.second.y, etak.second.z};
    metal.uroughness    = 0.2f;
    metal.vroughness    = 0.2f;
  } else {
    throw std::runtime_error("unknown pbrt bsdf filename " + fourier.bsdffile);
  }
}

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
static inline pair<vec3f, vec3f> parse_pbrt_subsurface(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> params = {
      // From "A Practical Model for Subsurface Light Transport"
      // Jensen, Marschner, Levoy, Hanrahan
      // Proc SIGGRAPH 2001
      {"Apple", {{2.29, 2.39, 1.97}, {0.0030, 0.0034, 0.046}}},
      {"Chicken1", {{0.15, 0.21, 0.38}, {0.015, 0.077, 0.19}}},
      {"Chicken2", {{0.19, 0.25, 0.32}, {0.018, 0.088, 0.20}}},
      {"Cream", {{7.38, 5.47, 3.15}, {0.0002, 0.0028, 0.0163}}},
      {"Ketchup", {{0.18, 0.07, 0.03}, {0.061, 0.97, 1.45}}},
      {"Marble", {{2.19, 2.62, 3.00}, {0.0021, 0.0041, 0.0071}}},
      {"Potato", {{0.68, 0.70, 0.55}, {0.0024, 0.0090, 0.12}}},
      {"Skimmilk", {{0.70, 1.22, 1.90}, {0.0014, 0.0025, 0.0142}}},
      {"Skin1", {{0.74, 0.88, 1.01}, {0.032, 0.17, 0.48}}},
      {"Skin2", {{1.09, 1.59, 1.79}, {0.013, 0.070, 0.145}}},
      {"Spectralon", {{11.6, 20.4, 14.9}, {0.00, 0.00, 0.00}}},
      {"Wholemilk", {{2.55, 3.21, 3.77}, {0.0011, 0.0024, 0.014}}},
      // From "Acquiring Scattering Properties of Participating Media by
      // Dilution",
      // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
      // Proc SIGGRAPH 2006
      {"Lowfat Milk", {{0.89187, 1.5136, 2.532}, {0.002875, 0.00575, 0.0115}}},
      {"Reduced Milk",
          {{2.4858, 3.1669, 4.5214}, {0.0025556, 0.0051111, 0.012778}}},
      {"Regular Milk",
          {{4.5513, 5.8294, 7.136}, {0.0015333, 0.0046, 0.019933}}},
      {"Espresso", {{0.72378, 0.84557, 1.0247}, {4.7984, 6.5751, 8.8493}}},
      {"Mint Mocha Coffee",
          {{0.31602, 0.38538, 0.48131}, {3.772, 5.8228, 7.82}}},
      {"Lowfat Soy Milk",
          {{0.30576, 0.34233, 0.61664}, {0.0014375, 0.0071875, 0.035937}}},
      {"Regular Soy Milk",
          {{0.59223, 0.73866, 1.4693}, {0.0019167, 0.0095833, 0.065167}}},
      {"Lowfat Chocolate Milk",
          {{0.64925, 0.83916, 1.1057}, {0.0115, 0.0368, 0.1564}}},
      {"Regular Chocolate Milk",
          {{1.4585, 2.1289, 2.9527}, {0.010063, 0.043125, 0.14375}}},
      {"Coke", {{8.9053e-05, 8.372e-05, 0}, {0.10014, 0.16503, 0.2468}}},
      {"Pepsi", {{6.1697e-05, 4.2564e-05, 0}, {0.091641, 0.14158, 0.20729}}},
      {"Sprite", {{6.0306e-06, 6.4139e-06, 6.5504e-06},
                     {0.001886, 0.0018308, 0.0020025}}},
      {"Gatorade",
          {{0.0024574, 0.003007, 0.0037325}, {0.024794, 0.019289, 0.008878}}},
      {"Chardonnay", {{1.7982e-05, 1.3758e-05, 1.2023e-05},
                         {0.010782, 0.011855, 0.023997}}},
      {"White Zinfandel", {{1.7501e-05, 1.9069e-05, 1.288e-05},
                              {0.012072, 0.016184, 0.019843}}},
      {"Merlot", {{2.1129e-05, 0, 0}, {0.11632, 0.25191, 0.29434}}},
      {"Budweiser Beer", {{2.4356e-05, 2.4079e-05, 1.0564e-05},
                             {0.011492, 0.024911, 0.057786}}},
      {"Coors Light Beer",
          {{5.0922e-05, 4.301e-05, 0}, {0.006164, 0.013984, 0.034983}}},
      {"Clorox",
          {{0.0024035, 0.0031373, 0.003991}, {0.0033542, 0.014892, 0.026297}}},
      {"Apple Juice",
          {{0.00013612, 0.00015836, 0.000227}, {0.012957, 0.023741, 0.052184}}},
      {"Cranberry Juice", {{0.00010402, 0.00011646, 7.8139e-05},
                              {0.039437, 0.094223, 0.12426}}},
      {"Grape Juice", {{5.382e-05, 0, 0}, {0.10404, 0.23958, 0.29325}}},
      {"Ruby Grapefruit Juice",
          {{0.011002, 0.010927, 0.011036}, {0.085867, 0.18314, 0.25262}}},
      {"White Grapefruit Juice",
          {{0.22826, 0.23998, 0.32748}, {0.0138, 0.018831, 0.056781}}},
      {"Shampoo",
          {{0.0007176, 0.0008303, 0.0009016}, {0.014107, 0.045693, 0.061717}}},
      {"Strawberry Shampoo",
          {{0.00015671, 0.00015947, 1.518e-05}, {0.01449, 0.05796, 0.075823}}},
      {"Head & Shoulders Shampoo",
          {{0.023805, 0.028804, 0.034306}, {0.084621, 0.15688, 0.20365}}},
      {"Lemon Tea Powder",
          {{0.040224, 0.045264, 0.051081}, {2.4288, 4.5757, 7.2127}}},
      {"Orange Powder", {{0.00015617, 0.00017482, 0.0001762},
                            {0.001449, 0.003441, 0.007863}}},
      {"Pink Lemonade Powder", {{0.00012103, 0.00013073, 0.00012528},
                                   {0.001165, 0.002366, 0.003195}}},
      {"Cappuccino Powder",
          {{1.8436, 2.5851, 2.1662}, {35.844, 49.547, 61.084}}},
      {"Salt Powder",
          {{0.027333, 0.032451, 0.031979}, {0.28415, 0.3257, 0.34148}}},
      {"Sugar Powder",
          {{0.00022272, 0.00025513, 0.000271}, {0.012638, 0.031051, 0.050124}}},
      {"Suisse Mocha Powder",
          {{2.7979, 3.5452, 4.3365}, {17.502, 27.004, 35.433}}},
      {"Pacific Ocean Surface Water", {{0.0001764, 0.00032095, 0.00019617},
                                          {0.031845, 0.031324, 0.030147}}},
  };
  return params.at(name);
}

// Get typename
static inline void parse_pbrt_typeparam(string_view& str, string& value) {
  auto saved = str;
  value      = "";
  auto pname = ""s, ptype = ""s;
  while (is_pbrt_param(str) && value == "") {
    parse_pbrt_nametype(str, pname, ptype);
    if (pname == "type") {
      parse_pbrt_param(str, ptype, value);
    } else {
      skip_pbrt_param(str);
    }
  }
  if (value == "") throw std::runtime_error("type not found");
  str = saved;
}

// Parse param and resolve constant textures
static inline void parse_pbrt_texture(string_view& str, const string& ptype,
    pbrt_textured3f&                              value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  parse_pbrt_param(str, ptype, value);
  if (value.texture == "") return;
  if (constant_values.find(value.texture) == constant_values.end()) return;
  value.value   = constant_values.at(value.texture);
  value.texture = "";
}
static inline void parse_pbrt_texture(string_view& str, const string& ptype,
    pbrt_textured1f&                              value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  parse_pbrt_param(str, ptype, value);
  if (value.texture == "") return;
  if (constant_values.find(value.texture) == constant_values.end()) return;
  auto col      = constant_values.at(value.texture);
  value.value   = (col.x + col.y + col.z) / 3;
  value.texture = "";
}

// Parse Material
static inline void parse_pbrt_material(string_view& str, const string& type,
    pbrt_material&                                value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  auto pname = ""s, ptype = ""s;
  if (type == "matte") {
    auto tvalue = pbrt_material::matte_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "sigma") {
        parse_pbrt_param(str, ptype, tvalue.sigma);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::matte;
    value.matte = tvalue;
  } else if (type == "mirror") {
    auto tvalue = pbrt_material::mirror_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_material::type_t::mirror;
    value.mirror = tvalue;
  } else if (type == "plastic") {
    auto tvalue = pbrt_material::plastic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_material::type_t::plastic;
    value.plastic = tvalue;
  } else if (type == "metal") {
    auto tvalue = pbrt_material::metal_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "k") {
        parse_pbrt_texture(str, ptype, tvalue.k, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::metal;
    value.metal = tvalue;
  } else if (type == "glass") {
    auto tvalue = pbrt_material::glass_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::glass;
    value.glass = tvalue;
  } else if (type == "translucent") {
    auto tvalue = pbrt_material::translucent_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "reflect") {
        parse_pbrt_texture(str, ptype, tvalue.reflect, constant_values);
      } else if (pname == "transmit") {
        parse_pbrt_texture(str, ptype, tvalue.transmit, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_material::type_t::translucent;
    value.translucent = tvalue;
  } else if (type == "uber") {
    auto tvalue = pbrt_material::uber_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "opacity") {
        parse_pbrt_texture(str, ptype, tvalue.opacity, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::uber;
    value.uber = tvalue;
  } else if (type == "disney") {
    auto tvalue = pbrt_material::disney_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "color") {
        parse_pbrt_texture(str, ptype, tvalue.color, constant_values);
      } else if (pname == "anisotropic") {
        parse_pbrt_texture(str, ptype, tvalue.anisotropic, constant_values);
      } else if (pname == "clearcoat") {
        parse_pbrt_texture(str, ptype, tvalue.clearcoat, constant_values);
      } else if (pname == "clearcoatgloss") {
        parse_pbrt_texture(str, ptype, tvalue.clearcoatgloss, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "metallic") {
        parse_pbrt_texture(str, ptype, tvalue.metallic, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "scatterdistance") {
        parse_pbrt_texture(str, ptype, tvalue.scatterdistance, constant_values);
      } else if (pname == "sheen") {
        parse_pbrt_texture(str, ptype, tvalue.sheen, constant_values);
      } else if (pname == "sheentint") {
        parse_pbrt_texture(str, ptype, tvalue.sheentint, constant_values);
      } else if (pname == "spectrans") {
        parse_pbrt_texture(str, ptype, tvalue.spectrans, constant_values);
      } else if (pname == "thin") {
        parse_pbrt_param(str, ptype, tvalue.thin);
      } else if (pname == "difftrans") {
        parse_pbrt_texture(str, ptype, tvalue.difftrans, constant_values);
      } else if (pname == "flatness") {
        parse_pbrt_texture(str, ptype, tvalue.flatness, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_material::type_t::disney;
    value.disney = tvalue;
  } else if (type == "hair") {
    auto tvalue = pbrt_material::hair_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "color") {
        parse_pbrt_texture(str, ptype, tvalue.color, constant_values);
      } else if (pname == "sigma_a") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_a, constant_values);
      } else if (pname == "eumelanin") {
        parse_pbrt_texture(str, ptype, tvalue.eumelanin, constant_values);
      } else if (pname == "pheomelanin") {
        parse_pbrt_texture(str, ptype, tvalue.pheomelanin, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "beta_m") {
        parse_pbrt_texture(str, ptype, tvalue.beta_m, constant_values);
      } else if (pname == "beta_n") {
        parse_pbrt_texture(str, ptype, tvalue.beta_n, constant_values);
      } else if (pname == "alpha") {
        parse_pbrt_texture(str, ptype, tvalue.alpha, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::hair;
    value.hair = tvalue;
  } else if (type == "kdsubsurface") {
    auto tvalue = pbrt_material::kdsubsurface_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "mfp") {
        parse_pbrt_texture(str, ptype, tvalue.mfp, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_material::type_t::kdsubsurface;
    value.kdsubsurface = tvalue;
  } else if (type == "mix") {
    auto tvalue = pbrt_material::mix_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "amount") {
        parse_pbrt_texture(str, ptype, tvalue.amount, constant_values);
      } else if (pname == "namedmaterial1") {
        parse_pbrt_param(str, ptype, tvalue.namedmaterial1);
      } else if (pname == "namedmaterial2") {
        parse_pbrt_param(str, ptype, tvalue.namedmaterial2);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::mix;
    value.mix  = tvalue;
  } else if (type == "fourier") {
    auto tvalue = pbrt_material::fourier_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "bsdffile") {
        parse_pbrt_param(str, ptype, tvalue.bsdffile);
        approximate_fourier_material(tvalue);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_material::type_t::fourier;
    value.fourier = tvalue;
  } else if (type == "substrate") {
    auto tvalue = pbrt_material::substrate_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type      = pbrt_material::type_t::substrate;
    value.substrate = tvalue;
  } else if (type == "subsurface") {
    auto tvalue = pbrt_material::subsurface_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "name") {
        parse_pbrt_param(str, ptype, tvalue.name);
        auto params    = parse_pbrt_subsurface(tvalue.name);
        tvalue.sigma_a = {params.second.x, params.second.y, params.second.z};
        tvalue.sigma_prime_s = {params.first.x, params.first.y, params.first.z};
      } else if (pname == "sigma_a") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_a, constant_values);
      } else if (pname == "sigma_prime_s") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_prime_s, constant_values);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_material::type_t::subsurface;
    value.subsurface = tvalue;
  } else {
    throw std::runtime_error("unknown Material " + type);
  }
}

// Parse Shape
static inline void parse_pbrt_shape(
    string_view& str, const string& type, pbrt_shape& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "trianglemesh") {
    auto tvalue = pbrt_shape::trianglemesh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "indices") {
        parse_pbrt_param(str, ptype, tvalue.indices);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "N") {
        parse_pbrt_param(str, ptype, tvalue.N);
      } else if (pname == "S") {
        parse_pbrt_param(str, ptype, tvalue.S);
      } else if (pname == "uv") {
        parse_pbrt_param(str, ptype, tvalue.uv);
      } else if (pname == "st") {
        parse_pbrt_param(str, ptype, tvalue.uv);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else if (pname == "shadowalpha") {
        parse_pbrt_param(str, ptype, tvalue.shadowalpha);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_shape::type_t::trianglemesh;
    value.trianglemesh = tvalue;
  } else if (type == "plymesh") {
    auto tvalue = pbrt_shape::plymesh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else if (pname == "shadowalpha") {
        parse_pbrt_param(str, ptype, tvalue.shadowalpha);
      } else if (pname == "discarddegenerateUVs") {
        // hack for some files
        auto value = false;
        parse_pbrt_param(str, ptype, value);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_shape::type_t::plymesh;
    value.plymesh = tvalue;
  } else if (type == "curve") {
    auto tvalue = pbrt_shape::curve_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "N") {
        parse_pbrt_param(str, ptype, tvalue.N);
      } else if (pname == "basis") {
        parse_pbrt_param(str, ptype, tvalue.basis);
      } else if (pname == "degree") {
        parse_pbrt_param(str, ptype, tvalue.degree);
      } else if (pname == "type") {
        parse_pbrt_param(str, ptype, tvalue.type);
      } else if (pname == "width") {
        auto width = 1.0f;
        parse_pbrt_param(str, ptype, width);
        tvalue.width0 = width;
        tvalue.width1 = width;
      } else if (pname == "width0") {
        parse_pbrt_param(str, ptype, tvalue.width0);
      } else if (pname == "width1") {
        parse_pbrt_param(str, ptype, tvalue.width1);
      } else if (pname == "splitdepth") {
        parse_pbrt_param(str, ptype, tvalue.splitdepth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_shape::type_t::curve;
    value.curve = tvalue;
  } else if (type == "loopsubdiv") {
    auto tvalue = pbrt_shape::loopsubdiv_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "indices") {
        parse_pbrt_param(str, ptype, tvalue.indices);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "levels") {
        parse_pbrt_param(str, ptype, tvalue.levels);
      } else if (pname == "nlevels") {
        parse_pbrt_param(str, ptype, tvalue.levels);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_shape::type_t::loopsubdiv;
    value.loopsubdiv = tvalue;
  } else if (type == "nurbs") {
    auto tvalue = pbrt_shape::nurbs_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "nu") {
        parse_pbrt_param(str, ptype, tvalue.nu);
      } else if (pname == "nv") {
        parse_pbrt_param(str, ptype, tvalue.nv);
      } else if (pname == "uknots") {
        parse_pbrt_param(str, ptype, tvalue.uknots);
      } else if (pname == "vknots") {
        parse_pbrt_param(str, ptype, tvalue.vknots);
      } else if (pname == "u0") {
        parse_pbrt_param(str, ptype, tvalue.u0);
      } else if (pname == "v0") {
        parse_pbrt_param(str, ptype, tvalue.v0);
      } else if (pname == "u1") {
        parse_pbrt_param(str, ptype, tvalue.u1);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "Pw") {
        parse_pbrt_param(str, ptype, tvalue.Pw);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_shape::type_t::nurbs;
    value.nurbs = tvalue;
  } else if (type == "sphere") {
    auto tvalue = pbrt_shape::sphere_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_shape::type_t::sphere;
    value.sphere = tvalue;
  } else if (type == "disk") {
    auto tvalue = pbrt_shape::disk_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "height") {
        parse_pbrt_param(str, ptype, tvalue.height);
      } else if (pname == "innerradius") {
        parse_pbrt_param(str, ptype, tvalue.innerradius);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_shape::type_t::disk;
    value.disk = tvalue;
  } else if (type == "cone") {
    auto tvalue = pbrt_shape::cone_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "height") {
        parse_pbrt_param(str, ptype, tvalue.height);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_shape::type_t::cone;
    value.cone = tvalue;
  } else if (type == "cylinder") {
    auto tvalue = pbrt_shape::cylinder_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_shape::type_t::cylinder;
    value.cylinder = tvalue;
  } else if (type == "hyperboloid") {
    auto tvalue = pbrt_shape::hyperboloid_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "p1") {
        parse_pbrt_param(str, ptype, tvalue.p1);
      } else if (pname == "p2") {
        parse_pbrt_param(str, ptype, tvalue.p2);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_shape::type_t::hyperboloid;
    value.hyperboloid = tvalue;
  } else if (type == "paraboloid") {
    auto tvalue = pbrt_shape::paraboloid_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_shape::type_t::paraboloid;
    value.paraboloid = tvalue;
  } else if (type == "heightfield") {
    auto tvalue = pbrt_shape::heightfield_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "nu") {
        parse_pbrt_param(str, ptype, tvalue.nu);
      } else if (pname == "nv") {
        parse_pbrt_param(str, ptype, tvalue.nv);
      } else if (pname == "Pz") {
        parse_pbrt_param(str, ptype, tvalue.Pz);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_shape::type_t::heightfield;
    value.heightfield = tvalue;
  } else {
    throw std::runtime_error("unknown Shape " + type);
  }
}

// Parse AreaLightSource
static inline void parse_pbrt_arealight(
    string_view& str, const string& type, pbrt_arealight& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "diffuse") {
    auto tvalue = pbrt_arealight::diffuse_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "twosided") {
        parse_pbrt_param(str, ptype, tvalue.twosided);
      } else if (pname == "samples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "nsamples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_arealight::type_t::diffuse;
    value.diffuse = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse LightSource
static inline void parse_pbrt_light(
    string_view& str, const string& type, pbrt_light& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "distant") {
    auto tvalue = pbrt_light::distant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else if (pname == "to") {
        parse_pbrt_param(str, ptype, tvalue.to);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_light::type_t::distant;
    value.distant = tvalue;
  } else if (type == "goniometric") {
    auto tvalue = pbrt_light::goniometric_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_light::type_t::goniometric;
    value.goniometric = tvalue;
  } else if (type == "infinite") {
    auto tvalue = pbrt_light::infinite_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "samples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "nsamples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_light::type_t::infinite;
    value.infinite = tvalue;
  } else if (type == "distant") {
    auto tvalue = pbrt_light::distant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_light::type_t::distant;
    value.distant = tvalue;
  } else if (type == "projection") {
    auto tvalue = pbrt_light::projection_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "fov") {
        parse_pbrt_param(str, ptype, tvalue.fov);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_light::type_t::projection;
    value.projection = tvalue;
  } else if (type == "spot") {
    auto tvalue = pbrt_light::spot_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else if (pname == "to") {
        parse_pbrt_param(str, ptype, tvalue.to);
      } else if (pname == "coneangle") {
        parse_pbrt_param(str, ptype, tvalue.coneangle);
      } else if (pname == "conedeltaangle") {
        parse_pbrt_param(str, ptype, tvalue.conedeltaangle);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_light::type_t::spot;
    value.spot = tvalue;
  } else if (type == "point") {
    auto tvalue = pbrt_light::point_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_light::type_t::point;
    value.point = tvalue;
  } else {
    throw std::runtime_error("unknown LightSource " + type);
  }
}

// Parse Medium
static inline void parse_pbrt_medium(
    string_view& str, const string& type, pbrt_medium& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "homogeneous") {
    auto tvalue = pbrt_medium::homogeneous_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "sigma_a") {
        parse_pbrt_param(str, ptype, tvalue.sigma_a);
      } else if (pname == "sigma_s") {
        parse_pbrt_param(str, ptype, tvalue.sigma_s);
      } else if (pname == "preset") {
        parse_pbrt_param(str, ptype, tvalue.preset);
      } else if (pname == "g") {
        parse_pbrt_param(str, ptype, tvalue.g);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_medium::type_t::homogeneous;
    value.homogeneous = tvalue;
  } else if (type == "heterogeneous") {
    auto tvalue = pbrt_medium::heterogeneous_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "sigma_a") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "sigma_s") {
        parse_pbrt_param(str, ptype, tvalue.sigma_s);
      } else if (pname == "preset") {
        parse_pbrt_param(str, ptype, tvalue.preset);
      } else if (pname == "g") {
        parse_pbrt_param(str, ptype, tvalue.g);
      } else if (pname == "p0") {
        parse_pbrt_param(str, ptype, tvalue.p0);
      } else if (pname == "p1") {
        parse_pbrt_param(str, ptype, tvalue.p1);
      } else if (pname == "nx") {
        parse_pbrt_param(str, ptype, tvalue.nx);
      } else if (pname == "ny") {
        parse_pbrt_param(str, ptype, tvalue.ny);
      } else if (pname == "nz") {
        parse_pbrt_param(str, ptype, tvalue.nz);
      } else if (pname == "density") {
        parse_pbrt_param(str, ptype, tvalue.density);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type          = pbrt_medium::type_t::heterogeneous;
    value.heterogeneous = tvalue;
  } else {
    throw std::runtime_error("unknown Medium " + type);
  }
}

// Load pbrt scene
void load_pbrt(const string& filename, pbrt_callbacks& cb, bool flipv) {
  // start laoding files
  auto  files = vector<file_holder>{};
  auto& file  = files.emplace_back();
  open_input_file(file, filename);

  // parsing stack
  auto stack    = vector<pbrt_context>{{}};
  auto object   = pbrt_object{};
  auto coordsys = unordered_map<string, pair<frame3f, frame3f>>{};

  // helpders
  auto set_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end = (frame3f)xform;
  };
  auto concat_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end *= (frame3f)xform;
  };

  // constant values
  unordered_map<string, pbrt_spectrum3f> constant_values = {};

  // parse command by command
  while (!files.empty()) {
    auto fs   = files.back().fs;
    auto line = ""s;
    auto cmd  = ""s;
    while (read_pbrt_line(fs, line)) {
      auto str = string_view{line};
      // get command
      parse_pbrt_command(str, cmd);
      if (cmd == "WorldBegin") {
        stack.push_back({});
      } else if (cmd == "WorldEnd") {
        stack.pop_back();
        if (stack.size() != 1) throw std::runtime_error("bad stack");
      } else if (cmd == "AttributeBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "AttributeEnd") {
        stack.pop_back();
      } else if (cmd == "TransformBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "TransformEnd") {
        stack.pop_back();
      } else if (cmd == "ObjectBegin") {
        parse_pbrt_value(str, object.name);
        stack.push_back(stack.back());
        cb.begin_object(object, stack.back());
      } else if (cmd == "ObjectEnd") {
        cb.end_object(object, stack.back());
        stack.pop_back();
        object = {};
      } else if (cmd == "ObjectInstance") {
        auto value = pbrt_object{};
        parse_pbrt_value(str, value.name);
        cb.object_instance(value, stack.back());
      } else if (cmd == "ActiveTransform") {
        auto value = ""s;
        parse_pbrt_command(str, value);
        if (value == "StartTime") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = false;
        } else if (value == "EndTime") {
          stack.back().active_transform_start = false;
          stack.back().active_transform_end   = true;
        } else if (value == "All") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = true;
        } else {
          throw std::runtime_error("bad active transform");
        }
      } else if (cmd == "Transform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        set_transform(stack.back(), xf);
      } else if (cmd == "ConcatTransform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        concat_transform(stack.back(), xf);
      } else if (cmd == "Scale") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), (mat4f)scaling_frame(v));
      } else if (cmd == "Translate") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), (mat4f)translation_frame(v));
      } else if (cmd == "Rotate") {
        auto v = zero4f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(),
            (mat4f)rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
      } else if (cmd == "LookAt") {
        auto from = zero3f, to = zero3f, up = zero3f;
        parse_pbrt_param(str, from);
        parse_pbrt_param(str, to);
        parse_pbrt_param(str, up);
        // from pbrt parser
        auto frame = lookat_frame(from, to, up, true);
        // frame.z = normalize(to-from);
        // frame.x = normalize(cross(frame.z,up));
        // frame.y = cross(frame.x,frame.z);
        // frame.o    = from;
        concat_transform(stack.back(), (mat4f)inverse(frame));
        stack.back().last_lookat_distance = length(from - to);
        // stack.back().focus = length(m.x - m.y);
      } else if (cmd == "ReverseOrientation") {
        stack.back().reverse = !stack.back().reverse;
      } else if (cmd == "CoordinateSystem") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        coordsys[name] = {
            stack.back().transform_start, stack.back().transform_end};
      } else if (cmd == "CoordSysTransform") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        if (coordsys.find(name) != coordsys.end()) {
          stack.back().transform_start = coordsys.at(name).first;
          stack.back().transform_end   = coordsys.at(name).second;
        }
      } else if (cmd == "Integrator") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_integrator{};
        parse_pbrt_integrator(str, type, value);
        cb.integrator(value, stack.back());
      } else if (cmd == "Sampler") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_sampler{};
        parse_pbrt_sampler(str, type, value);
        cb.sampler(value, stack.back());
      } else if (cmd == "PixelFilter") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_filter{};
        parse_pbrt_filter(str, type, value);
        cb.filter(value, stack.back());
      } else if (cmd == "Film") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_film{};
        parse_pbrt_film(str, type, value);
        cb.film(value, stack.back());
      } else if (cmd == "Accelerator") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_accelerator{};
        parse_pbrt_accelerator(str, type, value);
        cb.accelerator(value, stack.back());
      } else if (cmd == "Camera") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_camera{};
        parse_pbrt_camera(str, type, value);
        cb.camera(value, stack.back());
      } else if (cmd == "Texture") {
        auto name = ""s, comptype = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_value(str, comptype);
        parse_pbrt_value(str, type);
        auto value = pbrt_texture{};
        parse_pbrt_texture(str, type, value);
        if (type == "constant") {
          constant_values[name] = value.constant.value.value;
        }
        cb.texture(value, name, stack.back());
      } else if (cmd == "Material") {
        static auto material_id = 0;
        auto        type        = ""s;
        parse_pbrt_value(str, type);
        if (type == "") {
          stack.back().material = "";
        } else {
          auto value = pbrt_material{};
          auto name  = "unnamed_material_" + std::to_string(material_id++);
          parse_pbrt_material(str, type, value, constant_values);
          stack.back().material = name;
          cb.material(value, name, stack.back());
        }
      } else if (cmd == "MakeNamedMaterial") {
        auto name = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_typeparam(str, type);
        auto value = pbrt_material{};
        parse_pbrt_material(str, type, value, constant_values);
        cb.material(value, name, stack.back());
      } else if (cmd == "NamedMaterial") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        stack.back().material = name;
      } else if (cmd == "Shape") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_shape{};
        parse_pbrt_shape(str, type, value);
        cb.shape(value, stack.back());
      } else if (cmd == "AreaLightSource") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        static auto material_id = 0;
        auto        name = "unnamed_arealight_" + std::to_string(material_id++);
        auto        value = pbrt_arealight{};
        parse_pbrt_arealight(str, type, value);
        stack.back().arealight = name;
        cb.arealight(value, name, stack.back());
      } else if (cmd == "LightSource") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_light{};
        parse_pbrt_light(str, type, value);
        cb.light(value, stack.back());
      } else if (cmd == "MakeNamedMedium") {
        auto name = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_typeparam(str, type);
        auto value = pbrt_medium{};
        parse_pbrt_medium(str, type, value);
        cb.medium(value, name, stack.back());
      } else if (cmd == "MediumInterface") {
        auto interior = ""s, exterior = ""s;
        parse_pbrt_value(str, interior);
        parse_pbrt_value(str, exterior);
        stack.back().medium_interior = interior;
        stack.back().medium_exterior = exterior;
      } else if (cmd == "Include") {
        auto inputname = ""s;
        parse_pbrt_value(str, inputname);
        auto& file = files.emplace_back();
        open_input_file(file, fs::path(filename).parent_path() / inputname);
      } else {
        throw std::runtime_error("unknown command " + cmd);
      }
    }
    files.pop_back();
  }
}

// Load pbrt scene
bool read_pbrt_element(FILE* fs, pbrt_element& element, string& name,
    pbrt_element_data& data, vector<pbrt_context>& stack,
    pbrt_parser_state& state) {
  // helpders
  auto set_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end = (frame3f)xform;
  };
  auto concat_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end *= (frame3f)xform;
  };

  // init stack
  if (stack.empty()) stack.emplace_back();

  // parse command by command
  while (read_pbrt_line(fs, state.line)) {
    auto str = string_view{state.line};
    // get command
    auto cmd = ""s;
    parse_pbrt_command(str, cmd);
    if (cmd == "WorldBegin") {
      stack.push_back({});
    } else if (cmd == "WorldEnd") {
      stack.pop_back();
      if (stack.size() != 1) throw std::runtime_error("bad stack");
    } else if (cmd == "AttributeBegin") {
      stack.push_back(stack.back());
    } else if (cmd == "AttributeEnd") {
      stack.pop_back();
    } else if (cmd == "TransformBegin") {
      stack.push_back(stack.back());
    } else if (cmd == "TransformEnd") {
      stack.pop_back();
    } else if (cmd == "ObjectBegin") {
      parse_pbrt_value(str, state.object);
      stack.push_back(stack.back());
      element = pbrt_element::begin_object;
      name    = state.object;
      return true;
    } else if (cmd == "ObjectEnd") {
      element = pbrt_element::end_object;
      name    = state.object;
      stack.pop_back();
      state.object = {};
      return true;
    } else if (cmd == "ObjectInstance") {
      parse_pbrt_value(str, name);
      element = pbrt_element::object_instance;
      return true;
    } else if (cmd == "ActiveTransform") {
      auto value = ""s;
      parse_pbrt_command(str, value);
      if (value == "StartTime") {
        stack.back().active_transform_start = true;
        stack.back().active_transform_end   = false;
      } else if (value == "EndTime") {
        stack.back().active_transform_start = false;
        stack.back().active_transform_end   = true;
      } else if (value == "All") {
        stack.back().active_transform_start = true;
        stack.back().active_transform_end   = true;
      } else {
        throw std::runtime_error("bad active transform");
      }
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      set_transform(stack.back(), xf);
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      concat_transform(stack.back(), xf);
    } else if (cmd == "Scale") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(), (mat4f)scaling_frame(v));
    } else if (cmd == "Translate") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(), (mat4f)translation_frame(v));
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(),
          (mat4f)rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      parse_pbrt_param(str, from);
      parse_pbrt_param(str, to);
      parse_pbrt_param(str, up);
      // from pbrt parser
      auto frame = lookat_frame(from, to, up, true);
      // frame.z = normalize(to-from);
      // frame.x = normalize(cross(frame.z,up));
      // frame.y = cross(frame.x,frame.z);
      // frame.o    = from;
      concat_transform(stack.back(), (mat4f)inverse(frame));
      stack.back().last_lookat_distance = length(from - to);
      // stack.back().focus = length(m.x - m.y);
    } else if (cmd == "ReverseOrientation") {
      stack.back().reverse = !stack.back().reverse;
    } else if (cmd == "CoordinateSystem") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      state.coordsys[name] = {
          stack.back().transform_start, stack.back().transform_end};
    } else if (cmd == "CoordSysTransform") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      if (state.coordsys.find(name) != state.coordsys.end()) {
        stack.back().transform_start = state.coordsys.at(name).first;
        stack.back().transform_end   = state.coordsys.at(name).second;
      }
    } else if (cmd == "Integrator") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_integrator(str, type, data.intergrator);
      element = pbrt_element::integrator;
      return true;
    } else if (cmd == "Sampler") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_sampler(str, type, data.sampler);
      element = pbrt_element::sampler;
      return true;
    } else if (cmd == "PixelFilter") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_filter(str, type, data.filter);
      element = pbrt_element::filter;
      return true;
    } else if (cmd == "Film") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_film(str, type, data.film);
      element = pbrt_element::film;
      return true;
    } else if (cmd == "Accelerator") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_accelerator(str, type, data.accelerator);
      element = pbrt_element::accelerator;
      return true;
    } else if (cmd == "Camera") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_camera(str, type, data.camera);
      element = pbrt_element::camera;
      return true;
    } else if (cmd == "Texture") {
      auto comptype = ""s, type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_value(str, comptype);
      parse_pbrt_value(str, type);
      parse_pbrt_texture(str, type, data.texture);
      if (type == "constant") {
        state.constant_values[name] = data.texture.constant.value.value;
      }
      element = pbrt_element::texture;
      return true;
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        type        = ""s;
      parse_pbrt_value(str, type);
      if (type == "") {
        stack.back().material = "";
      } else {
        name = "unnamed_material_" + std::to_string(material_id++);
        parse_pbrt_material(str, type, data.material, state.constant_values);
        stack.back().material = name;
        element               = pbrt_element::material;
        return true;
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_typeparam(str, type);
      parse_pbrt_material(str, type, data.material, state.constant_values);
      element = pbrt_element::material;
      return true;
    } else if (cmd == "NamedMaterial") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      stack.back().material = name;
    } else if (cmd == "Shape") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_shape(str, type, data.shape);
      element = pbrt_element::shape;
      return true;
    } else if (cmd == "AreaLightSource") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      static auto material_id = 0;
      name = "unnamed_arealight_" + std::to_string(material_id++);
      parse_pbrt_arealight(str, type, data.arealight);
      stack.back().arealight = name;
      element                = pbrt_element::arealight;
      return true;
    } else if (cmd == "LightSource") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_light(str, type, data.light);
      element = pbrt_element::light;
      return true;
    } else if (cmd == "MakeNamedMedium") {
      auto type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_typeparam(str, type);
      parse_pbrt_medium(str, type, data.medium);
      element = pbrt_element::medium;
      return true;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      parse_pbrt_value(str, interior);
      parse_pbrt_value(str, exterior);
      stack.back().medium_interior = interior;
      stack.back().medium_exterior = exterior;
    } else if (cmd == "Include") {
      parse_pbrt_value(str, name);
      element = pbrt_element::include;
      return true;
    } else {
      throw std::runtime_error("unknown command " + cmd);
    }
  }
  return false;
}

}  // namespace yocto
