//
// Implementation for Yocto/PbrtIO.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include "yocto_pbrtio.h"

#include <fast_float/fast_float.h>

#include "yocto_modelio.h"

#define _USE_MATH_DEFINES
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string_view;
using std::unordered_map;
using namespace std::string_literals;
using namespace std::string_view_literals;
using byte = unsigned char;
using std::cos;
using std::sin;

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8    = std::filesystem::u8path(filename);
  auto str_mode = string{mode};
  auto wmode    = std::wstring(str_mode.begin(), str_mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Load a text file
static bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = "cannot open " + filename;
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    error = "cannot read " + filename;
    return false;
  }
  fclose(fs);
  return true;
}

// Save a text file
static bool save_text(
    const string& filename, const string& str, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (!fs) {
    error = "cannot create " + filename;
    return false;
  }
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    error = "cannot write " + filename;
    return false;
  }
  fclose(fs);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
static std::filesystem::path make_path(const string& filename) {
  return std::filesystem::u8path(filename);
}

// Get directory name (not including /)
static string path_dirname(const string& filename) {
  return make_path(filename).parent_path().generic_u8string();
}

// Get filename without directory.
static string path_filename(const string& filename) {
  return make_path(filename).filename().generic_u8string();
}

// Joins paths
static string path_join(const string& patha, const string& pathb) {
  return (make_path(patha) / make_path(pathb)).generic_u8string();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ARRAY MATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

static array<float, 3> neg(const array<float, 3>& a) {
  return {-a[0], -a[1], -a[2]};
}
static array<float, 3> add(const array<float, 3>& a, const array<float, 3>& b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
static array<float, 3> sub(const array<float, 3>& a, const array<float, 3>& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
static array<float, 3> mul(const array<float, 3>& a, float b) {
  return {a[0] * b, a[1] * b, a[2] * b};
}
static array<float, 3> div(const array<float, 3>& a, float b) {
  return {a[0] / b, a[1] / b, a[2] / b};
}

static float dot(const array<float, 3>& a, const array<float, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static array<float, 3> cross(
    const array<float, 3>& a, const array<float, 3>& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]};
}
static float length(const array<float, 3>& a) { return sqrt(dot(a, a)); }
static array<float, 3> normalize(const array<float, 3>& a) {
  auto l = length(a);
  return (l != 0) ? div(a, l) : a;
}

static array<array<float, 3>, 4> lookat_frame(const array<float, 3>& eye,
    const array<float, 3>& center, const array<float, 3>& up,
    bool inv_xz = false) {
  auto w = normalize(sub(eye, center));
  auto u = normalize(cross(up, w));
  auto v = normalize(cross(w, u));
  if (inv_xz) {
    w = neg(w);
    u = neg(u);
  }
  return {u, v, w, eye};
}

static array<array<float, 3>, 3> mul(
    const array<array<float, 3>, 3>& a, float b) {
  return {mul(a[0], b), mul(a[1], b), mul(a[2], b)};
}
static array<float, 3> mul(
    const array<array<float, 3>, 3>& a, const array<float, 3>& b) {
  return add(mul(a[0], b[0]), add(mul(a[1], b[1]), mul(a[2], b[2])));
}
static array<array<float, 3>, 4> mul(
    const array<array<float, 3>, 4>& a, const array<array<float, 3>, 4>& b) {
  auto al = array<array<float, 3>, 3>{a[0], a[1], a[2]};
  return {
      mul(al, b[0]), mul(al, b[1]), mul(al, b[2]), add(mul(al, b[3]), a[3])};
}

static array<array<float, 3>, 4> translation_frame(const array<float, 3>& a) {
  return {array<float, 3>{1, 0, 0}, array<float, 3>{0, 1, 0},
      array<float, 3>{0, 0, 1}, a};
}
static array<array<float, 3>, 4> scaling_frame(const array<float, 3>& a) {
  return {array<float, 3>{a[0], 0, 0}, array<float, 3>{0, a[1], 0},
      array<float, 3>{0, 0, a[2]}, array<float, 3>{0, 0, 0}};
}
static array<array<float, 3>, 4> rotation_frame(
    const array<float, 3>& axis, float angle) {
  auto s = std::sin(angle), c = std::cos(angle);
  auto vv = normalize(axis);
  return {array<float, 3>{c + (1 - c) * vv[0] * vv[0],
              (1 - c) * vv[0] * vv[1] + s * vv[2],
              (1 - c) * vv[0] * vv[2] - s * vv[1]},
      array<float, 3>{(1 - c) * vv[0] * vv[1] - s * vv[2],
          c + (1 - c) * vv[1] * vv[1], (1 - c) * vv[1] * vv[2] + s * vv[0]},
      array<float, 3>{(1 - c) * vv[0] * vv[2] + s * vv[1],
          (1 - c) * vv[1] * vv[2] - s * vv[0], c + (1 - c) * vv[2] * vv[2]},
      array<float, 3>{0, 0, 0}};
}

static array<array<float, 3>, 3> transpose(const array<array<float, 3>, 3>& a) {
  return {
      array<float, 3>{a[0][0], a[1][0], a[2][0]},
      array<float, 3>{a[0][1], a[1][1], a[2][1]},
      array<float, 3>{a[0][2], a[1][2], a[2][2]},
  };
}
static float determinant(const array<array<float, 3>, 3>& a) {
  return dot(a[0], cross(a[1], a[2]));
}
static array<array<float, 3>, 3> adjoint(const array<array<float, 3>, 3>& a) {
  return transpose(array<array<float, 3>, 3>{
      cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])});
}
static array<array<float, 3>, 3> inverse(const array<array<float, 3>, 3>& a) {
  return mul(adjoint(a), (1 / determinant(a)));
}

static array<array<float, 3>, 4> inverse(
    const array<array<float, 3>, 4>& a, bool non_rigid = false) {
  auto m    = array<array<float, 3>, 3>{a[0], a[1], a[2]};
  auto minv = non_rigid ? inverse(m) : transpose(m);
  return {minv[0], minv[1], minv[2], neg(mul(minv, a[3]))};
}

// frame/mat conversion
static array<array<float, 3>, 4> mat_to_frame(
    const array<array<float, 4>, 4>& m) {
  return {array<float, 3>{m[0][0], m[0][1], m[0][2]},
      array<float, 3>{m[1][0], m[1][1], m[1][2]},
      array<float, 3>{m[2][0], m[2][1], m[2][2]},
      array<float, 3>{m[3][0], m[3][1], m[3][2]}};
}
static array<array<float, 4>, 4> frame_to_mat(
    const array<array<float, 3>, 4>& f) {
  return {array<float, 4>{f[0][0], f[0][1], f[0][2], 0},
      array<float, 4>{f[1][0], f[1][1], f[1][2], 0},
      array<float, 4>{f[2][0], f[2][1], f[2][2], 0},
      array<float, 4>{f[3][0], f[3][1], f[3][2], 1}};
}

static array<float, 12> flatten(const array<array<float, 3>, 4>& a) {
  return (const array<float, 12>&)a;
}
static array<array<float, 3>, 4> unflatten(const array<float, 12>& a) {
  return (const array<array<float, 3>, 4>&)a;
}

static array<array<float, 4>, 4> unflatten(const array<float, 16>& a) {
  return (const array<array<float, 4>, 4>&)a;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR STL
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <typename T, size_t N>
struct hash<array<T, N>> {
  size_t operator()(const array<T, N>& value) const {
    const std::hash<float> hasher = std::hash<T>();
    auto                   h      = (size_t)0;
    for (auto item : value) {
      h ^= hasher(item) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
static void format_value(string& str, const string& value) { str += value; }
static void format_value(string& str, const char* value) { str += value; }
template <typename T>
static void format_value(string& str, T value) {
  auto buffer = array<char, 64>{};
  if constexpr (std::is_same_v<T, float>) {
#ifdef _WIN32
    auto result = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value);
    str.append(buffer.data(), result.ptr);
#else
    auto len = snprintf(buffer.data(), buffer.size(), "%.9g", value);
    str.append(buffer.data(), buffer.data() + len);
#endif
  } else if constexpr (std::is_same_v<T, double>) {
#ifdef _WIN32
    auto result = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value);
    str.append(buffer.data(), result.ptr);
#else
    auto len = snprintf(buffer.data(), buffer.size(), "%.17g", value);
    str.append(buffer.data(), buffer.data() + len);
#endif
  } else {
    auto result = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value);
    str.append(buffer.data(), result.ptr);
  }
}
template <typename T, size_t N>
static void format_value(string& str, const array<T, N>& value) {
  for (auto i = 0; i < N; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}

// Foramt to file
static void format_values(string& str, string_view fmt) {
  auto pos = fmt.find("{}"sv);
  if (pos != string::npos) throw std::invalid_argument("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static void format_values(
    string& str, string_view fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}"sv);
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

static bool is_newline(char c) { return c == '\r' || c == '\n'; }
static bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

static void remove_comment(
    string_view& str, char comment_char = '#', bool handle_quotes = false) {
  if (!handle_quotes) {
    while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
    str.remove_suffix(cpy.size());
  } else {
    while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
    auto cpy       = str;
    auto in_string = false;
    while (!cpy.empty()) {
      if (cpy.front() == '"') in_string = !in_string;
      if (cpy.front() == comment_char && !in_string) break;
      cpy.remove_prefix(1);
    }
    str.remove_suffix(cpy.size());
  }
}

// Read a line
static bool read_line(string_view& str, string_view& line) {
  if (str.empty()) return false;
  auto data = str.data();
  auto size = (size_t)0;
  while (!str.empty()) {
    if (str.front() == '\n') {
      str.remove_prefix(1);
      size++;
      break;
    } else {
      str.remove_prefix(1);
      size++;
    }
  }
  line = {data, size};
  return true;
}

// Parse values from a string
[[nodiscard]] static bool parse_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) return false;
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
    if (str.front() != '"') return false;
    str.remove_prefix(1);
    if (str.empty()) return false;
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
    if (cpy.empty()) return false;
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    str.remove_prefix(1);
  }
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
template <typename T>
[[nodiscard]] static bool parse_value(string_view& str, T& value) {
  skip_whitespace(str);
  if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
    auto result = fast_float::from_chars(
        str.data(), str.data() + str.size(), value);
    if (result.ptr == str.data()) return false;
    str.remove_prefix(result.ptr - str.data());
  } else {
    auto result = std::from_chars(str.data(), str.data() + str.size(), value);
    if (result.ptr == str.data()) return false;
    str.remove_prefix(result.ptr - str.data());
  }
  return true;
}
template <typename T, size_t N>
[[nodiscard]] static bool parse_value(string_view& str, array<T, N>& value) {
  for (auto i = 0; i < N; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

template <typename T>
[[nodiscard]] static bool read_value(string_view& str, T& value) {
  if (str.size() < sizeof(value)) return false;
  memcpy(&value, str.data(), sizeof(T));
  str.remove_prefix(sizeof(T));
  return true;
}

// Write data from a file
template <typename T>
static void write_value(vector<byte>& data, const T& value) {
  if constexpr (sizeof(T) == 1) {
    data.push_back(*(byte*)&value);
  } else {
    data.insert(data.end(), (byte*)(&value), (byte*)(&value) + sizeof(T));
  }
}

template <typename T>
static T swap_endian(T value) {
  // https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  static_assert(sizeof(char) == 1, "sizeof(char) == 1");
  union {
    T             value;
    unsigned char bytes[sizeof(T)];
  } source, dest;
  source.value = value;
  for (auto k = (size_t)0; k < sizeof(T); k++)
    dest.bytes[k] = source.bytes[sizeof(T) - k - 1];
  return dest.value;
}

template <typename T>
[[nodiscard]] static bool read_value(
    string_view& fs, T& value, bool big_endian) {
  if (!read_value(fs, value)) return false;
  if (big_endian) value = swap_endian(value);
  return true;
}

template <typename T>
static void write_value(vector<byte>& data, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  return write_value(data, value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt type
enum struct pbrt_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color,
  point2, vector2, spectrum
  // clang-format on
};

// Pbrt value
struct pbrt_value {
  string                  name     = "";
  pbrt_type               type     = pbrt_type::real;
  int                     value1i  = 0;
  float                   value1f  = 0;
  array<float, 2>         value2f  = {0, 0};
  array<float, 3>         value3f  = {0, 0, 0};
  bool                    value1b  = false;
  string                  value1s  = "";
  vector<float>           vector1f = {};
  vector<array<float, 2>> vector2f = {};
  vector<array<float, 3>> vector3f = {};
  vector<int>             vector1i = {};
};

// Pbrt command
struct pbrt_command {
  string                    name   = "";
  string                    type   = "";
  vector<pbrt_value>        values = {};
  array<array<float, 3>, 4> frame  = {array<float, 3>{1, 0, 0},
      array<float, 3>{0, 1, 0}, array<float, 3>{0, 0, 1},
      array<float, 3>{0, 0, 0}};
  array<array<float, 3>, 4> frend  = {array<float, 3>{1, 0, 0},
      array<float, 3>{0, 1, 0}, array<float, 3>{0, 0, 1},
      array<float, 3>{0, 0, 0}};
};

// get pbrt value
static bool get_pbrt_value(const pbrt_value& pbrt, string& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val = pbrt.value1s;
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(const pbrt_value& pbrt, bool& val) {
  if (pbrt.type == pbrt_type::boolean) {
    val = pbrt.value1b;
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(const pbrt_value& pbrt, int& val) {
  if (pbrt.type == pbrt_type::integer) {
    val = pbrt.value1i;
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(const pbrt_value& pbrt, float& val) {
  if (pbrt.type == pbrt_type::real) {
    val = pbrt.value1f;
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(const pbrt_value& pbrt, array<float, 3>& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    val = pbrt.value3f;
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    val = array<float, 3>{pbrt.value1f, pbrt.value1f, pbrt.value1f};
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(
    const pbrt_value& pbrt, vector<array<float, 2>>& val) {
  if (pbrt.type == pbrt_type::point2 || pbrt.type == pbrt_type::vector2) {
    if (!pbrt.vector2f.empty()) {
      val = pbrt.vector2f;
    } else {
      val = {pbrt.value2f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 2) != 0) return false;
    val.resize(pbrt.vector1f.size() / 2);
    for (auto i = 0; i < (int)val.size(); i++)
      val[i] = {pbrt.vector1f[i * 2 + 0], pbrt.vector1f[i * 2 + 1]};
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(
    const pbrt_value& pbrt, vector<array<float, 3>>& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    if (!pbrt.vector3f.empty()) {
      val = pbrt.vector3f;
    } else {
      val = {pbrt.value3f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 3) != 0) return false;
    val.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < (int)val.size(); i++)
      val[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}

static bool get_pbrt_value(const pbrt_value& pbrt, vector<array<int, 3>>& val) {
  if (pbrt.type == pbrt_type::integer) {
    if (pbrt.vector1i.empty() || (pbrt.vector1i.size() % 3) != 0) return false;
    val.resize(pbrt.vector1i.size() / 3);
    for (auto i = 0; i < (int)val.size(); i++)
      val[i] = {pbrt.vector1i[i * 3 + 0], pbrt.vector1i[i * 3 + 1],
          pbrt.vector1i[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}
static bool get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = 0;
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
static bool get_pbrt_value(
    const pbrt_value& pbrt, pair<array<float, 3>, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = {0, 0, 0};
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
template <typename T>
static bool get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& val) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_pbrt_value(p, val);
    }
  }
  return true;
}

// pbrt value construction
static pbrt_value make_pbrt_value(
    const string& name, const string& val, pbrt_type type = pbrt_type::string) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(
    const string& name, bool val, pbrt_type type = pbrt_type::boolean) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(
    const string& name, int val, pbrt_type type = pbrt_type::integer) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(
    const string& name, float val, pbrt_type type = pbrt_type::real) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(const string& name,
    const array<float, 3>& val, pbrt_type type = pbrt_type::color) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(const string& name,
    const vector<array<float, 2>>& val, pbrt_type type = pbrt_type::point2) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(const string& name,
    const vector<array<float, 3>>& val, pbrt_type type = pbrt_type::point) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = val;
  return pbrt;
}
static pbrt_value make_pbrt_value(const string& name,
    const vector<array<int, 3>>& val, pbrt_type type = pbrt_type::integer) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {&val.front()[0], &val.front()[0] + val.size() * 3};
  return pbrt;
}

[[maybe_unused]] static void remove_pbrt_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
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
static bool read_pbrt_cmdline(string_view& str, string& cmd) {
  cmd.clear();
  auto found = false;
  auto copy  = str;
  auto line  = string_view{};
  while (read_line(str, line)) {
    // line
    remove_comment(line, '#', true);
    skip_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        str = copy;
        return true;
      } else {
        found = true;
      }
    } else if (!found) {
      return false;
    }
    cmd += line;
    cmd += " ";
    copy = str;
  }
  return found;
}

// parse a quoted string
[[nodiscard]] static bool parse_command(string_view& str, string& value) {
  skip_whitespace(str);
  if (!isalpha((int)str.front())) return false;
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == string_view::npos) {
    value.assign(str);
    str.remove_prefix(str.size());
  } else {
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
  }
  return true;
}

// parse pbrt value with optional parens
template <typename T>
[[nodiscard]] static bool parse_param(string_view& str, T& value) {
  skip_whitespace(str);
  auto parens = !str.empty() && str.front() == '[';
  if (parens) str.remove_prefix(1);
  if (!parse_value(str, value)) return false;
  if (!str.data()) return false;
  if (parens) {
    skip_whitespace(str);
    if (!str.empty() && str.front() == '[') return false;
    str.remove_prefix(1);
  }
  return true;
}

// parse a quoted string
[[nodiscard]] static bool parse_nametype(
    string_view& str_, string& name, string& type) {
  auto value = ""s;
  if (!parse_value(str_, value)) return false;
  if (str_.empty()) return false;
  auto str  = string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == string_view::npos) return false;
  type = string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == string_view::npos) return false;
  str.remove_prefix(pos2);
  name = string(str);
  return true;
}

static pair<array<float, 3>, array<float, 3>> get_etak(const string& name) {
  static const unordered_map<string, pair<array<float, 3>, array<float, 3>>>
      metal_ior_table = {
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

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
[[maybe_unused]] static pair<array<float, 3>, array<float, 3>> get_subsurface(
    const string& name) {
  static const unordered_map<string, pair<array<float, 3>, array<float, 3>>>
      params = {
          // From "A Practical Model for Subsurface Light Transport"
          // Jensen, Marschner, Levoy, Hanrahan
          // Proc SIGGRAPH 2001
          {"Apple", {{2.29f, 2.39f, 1.97f}, {0.0030f, 0.0034f, 0.046f}}},
          {"Chicken1", {{0.15f, 0.21f, 0.38f}, {0.015f, 0.077f, 0.19f}}},
          {"Chicken2", {{0.19f, 0.25f, 0.32f}, {0.018f, 0.088f, 0.20f}}},
          {"Cream", {{7.38f, 5.47f, 3.15f}, {0.0002f, 0.0028f, 0.0163f}}},
          {"Ketchup", {{0.18f, 0.07f, 0.03f}, {0.061f, 0.97f, 1.45f}}},
          {"Marble", {{2.19f, 2.62f, 3.00f}, {0.0021f, 0.0041f, 0.0071f}}},
          {"Potato", {{0.68f, 0.70f, 0.55f}, {0.0024f, 0.0090f, 0.12f}}},
          {"Skimmilk", {{0.70f, 1.22f, 1.90f}, {0.0014f, 0.0025f, 0.0142f}}},
          {"Skin1", {{0.74f, 0.88f, 1.01f}, {0.032f, 0.17f, 0.48f}}},
          {"Skin2", {{1.09f, 1.59f, 1.79f}, {0.013f, 0.070f, 0.145f}}},
          {"Spectralon", {{11.6f, 20.4f, 14.9f}, {0.00f, 0.00f, 0.00f}}},
          {"Wholemilk", {{2.55f, 3.21f, 3.77f}, {0.0011f, 0.0024f, 0.014f}}},
          // From "Acquiring Scattering Properties of Participating Media by
          // Dilution",
          // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
          // Proc SIGGRAPH 2006
          {"Lowfat Milk",
              {{0.89187f, 1.5136f, 2.532f}, {0.002875f, 0.00575f, 0.0115f}}},
          {"Reduced Milk", {{2.4858f, 3.1669f, 4.5214f},
                               {0.0025556f, 0.0051111f, 0.012778f}}},
          {"Regular Milk",
              {{4.5513f, 5.8294f, 7.136f}, {0.0015333f, 0.0046f, 0.019933f}}},
          {"Espresso",
              {{0.72378f, 0.84557f, 1.0247f}, {4.7984f, 6.5751f, 8.8493f}}},
          {"Mint Mocha Coffee",
              {{0.31602f, 0.38538f, 0.48131f}, {3.772f, 5.8228f, 7.82f}}},
          {"Lowfat Soy Milk", {{0.30576f, 0.34233f, 0.61664f},
                                  {0.0014375f, 0.0071875f, 0.035937f}}},
          {"Regular Soy Milk", {{0.59223f, 0.73866f, 1.4693f},
                                   {0.0019167f, 0.0095833f, 0.065167f}}},
          {"Lowfat Chocolate Milk",
              {{0.64925f, 0.83916f, 1.1057f}, {0.0115f, 0.0368f, 0.1564f}}},
          {"Regular Chocolate Milk",
              {{1.4585f, 2.1289f, 2.9527f}, {0.010063f, 0.043125f, 0.14375f}}},
          {"Coke",
              {{8.9053e-05f, 8.372e-05f, 0.0f}, {0.10014f, 0.16503f, 0.2468f}}},
          {"Pepsi", {{6.1697e-05f, 4.2564e-05f, 0.0f},
                        {0.091641f, 0.14158f, 0.20729f}}},
          {"Sprite", {{6.0306e-06f, 6.4139e-06f, 6.5504e-06f},
                         {0.001886f, 0.0018308f, 0.0020025f}}},
          {"Gatorade", {{0.0024574f, 0.003007f, 0.0037325f},
                           {0.024794f, 0.019289f, 0.008878f}}},
          {"Chardonnay", {{1.7982e-05f, 1.3758e-05f, 1.2023e-05f},
                             {0.010782f, 0.011855f, 0.023997f}}},
          {"White Zinfandel", {{1.7501e-05f, 1.9069e-05f, 1.288e-05f},
                                  {0.012072f, 0.016184f, 0.019843f}}},
          {"Merlot",
              {{2.1129e-05f, 0.0f, 0.0f}, {0.11632f, 0.25191f, 0.29434f}}},
          {"Budweiser Beer", {{2.4356e-05f, 2.4079e-05f, 1.0564e-05f},
                                 {0.011492f, 0.024911f, 0.057786f}}},
          {"Coors Light Beer", {{5.0922e-05f, 4.301e-05f, 0.0f},
                                   {0.006164f, 0.013984f, 0.034983f}}},
          {"Clorox", {{0.0024035f, 0.0031373f, 0.003991f},
                         {0.0033542f, 0.014892f, 0.026297f}}},
          {"Apple Juice", {{0.00013612f, 0.00015836f, 0.000227f},
                              {0.012957f, 0.023741f, 0.052184f}}},
          {"Cranberry Juice", {{0.00010402f, 0.00011646f, 7.8139e-05f},
                                  {0.039437f, 0.094223f, 0.12426f}}},
          {"Grape Juice",
              {{5.382e-05f, 0.0f, 0.0f}, {0.10404f, 0.23958f, 0.29325f}}},
          {"Ruby Grapefruit Juice", {{0.011002f, 0.010927f, 0.011036f},
                                        {0.085867f, 0.18314f, 0.25262f}}},
          {"White Grapefruit Juice", {{0.22826f, 0.23998f, 0.32748f},
                                         {0.0138f, 0.018831f, 0.056781f}}},
          {"Shampoo", {{0.0007176f, 0.0008303f, 0.0009016f},
                          {0.014107f, 0.045693f, 0.061717f}}},
          {"Strawberry Shampoo", {{0.00015671f, 0.00015947f, 1.518e-05f},
                                     {0.01449f, 0.05796f, 0.075823f}}},
          {"Head & Shoulders Shampoo", {{0.023805f, 0.028804f, 0.034306f},
                                           {0.084621f, 0.15688f, 0.20365f}}},
          {"Lemon Tea Powder",
              {{0.040224f, 0.045264f, 0.051081f}, {2.4288f, 4.5757f, 7.2127f}}},
          {"Orange Powder", {{0.00015617f, 0.00017482f, 0.0001762f},
                                {0.001449f, 0.003441f, 0.007863f}}},
          {"Pink Lemonade Powder", {{0.00012103f, 0.00013073f, 0.00012528f},
                                       {0.001165f, 0.002366f, 0.003195f}}},
          {"Cappuccino Powder",
              {{1.8436f, 2.5851f, 2.1662f}, {35.844f, 49.547f, 61.084f}}},
          {"Salt Powder", {{0.027333f, 0.032451f, 0.031979f},
                              {0.28415f, 0.3257f, 0.34148f}}},
          {"Sugar Powder", {{0.00022272f, 0.00025513f, 0.000271f},
                               {0.012638f, 0.031051f, 0.050124f}}},
          {"Suisse Mocha Powder",
              {{2.7979f, 3.5452f, 4.3365f}, {17.502f, 27.004f, 35.433f}}},
          {"Pacific Ocean Surface Water",
              {{0.0001764f, 0.00032095f, 0.00019617f},
                  {0.031845f, 0.031324f, 0.030147f}}},
      };
  return params.at(name);
}

template <typename T, typename V>
[[nodiscard]] static bool parse_pvalues(string_view& str, T& value, V& values) {
  values.clear();
  skip_whitespace(str);
  if (str.empty()) return false;
  if (str.front() == '[') {
    str.remove_prefix(1);
    skip_whitespace(str);
    if (str.empty()) return false;
    while (!str.empty()) {
      auto& val = values.empty() ? value : values.emplace_back();
      if (!parse_value(str, val)) return false;
      if (!str.data()) return false;
      skip_whitespace(str);
      if (str.empty()) break;
      if (str.front() == ']') break;
      if (values.empty()) values.push_back(value);
    }
    if (str.empty()) return false;
    if (str.front() != ']') return false;
    str.remove_prefix(1);
  } else {
    if (!parse_value(str, value)) return false;
  }
  return true;
}

// Approximate color of blackbody radiation from wavelength in nm.
static array<float, 3> blackbody_to_rgb_(float temperature) {
  // clamp to valid range
  auto t = std::min(std::max(temperature, 1667.0f), 25000.0f) / 1000.0f;
  // compute x
  auto x = 0.0f;
  if (temperature < 4000.0f) {
    x = -0.2661239f * 1 / (t * t * t) - 0.2343589f * 1 / (t * t) +
        0.8776956f * (1 / t) + 0.179910f;
  } else {
    x = -3.0258469f * 1 / (t * t * t) + 2.1070379f * 1 / (t * t) +
        0.2226347f * (1 / t) + 0.240390f;
  }
  // compute y
  auto y = 0.0f;
  if (temperature < 2222.0f) {
    y = -1.1063814f * (x * x * x) - 1.34811020f * (x * x) + 2.18555832f * x -
        0.20219683f;
  } else if (temperature < 4000.0f) {
    y = -0.9549476f * (x * x * x) - 1.37418593f * (x * x) + 2.09137015f * x -
        0.16748867f;
  } else {
    y = +3.0817580f * (x * x * x) - 5.87338670f * (x * x) + 3.75112997f * x -
        0.37001483f;
  }
  auto xyY = array<float, 3>{x, y, 1};
  auto xyz = array<float, 3>{xyY[0] * xyY[2] / xyY[1], xyY[2],
      (1 - xyY[0] - xyY[1]) * xyY[2] / xyY[1]};
  auto rgb = add(mul(array<float, 3>{+3.2406f, -0.9689f, +0.0557f}, xyz[0]),
      add(mul(array<float, 3>{-1.5372f, +1.8758f, -0.2040f}, xyz[1]),
          mul(array<float, 3>{-0.4986f, +0.0415f, +1.0570f}, xyz[2])));
  return rgb;
}

[[nodiscard]] static bool parse_params(
    string_view& str, vector<pbrt_value>& values) {
  auto starts_with = [](string_view value, string_view prefix) {
    if (prefix.size() > value.size()) return false;
    return value.rfind(prefix, 0) == 0;
  };
  auto ends_with = [](string_view value, string_view postfix) {
    if (postfix.size() > value.size()) return false;
    return std::equal(postfix.rbegin(), postfix.rend(), value.rbegin());
  };

  values.clear();
  skip_whitespace(str);
  while (!str.empty()) {
    auto& value = values.emplace_back();
    auto  type  = ""s;
    if (!parse_nametype(str, value.name, type)) return false;
    skip_whitespace(str);
    if (str.empty()) return false;
    if (type == "float") {
      value.type = pbrt_type::real;
      if (!parse_pvalues(str, value.value1f, value.vector1f)) return false;
    } else if (type == "integer") {
      value.type = pbrt_type::integer;
      if (!parse_pvalues(str, value.value1i, value.vector1i)) return false;
    } else if (type == "string") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::string;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::boolean;
      if (!parse_pvalues(str, value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::texture;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "point" || type == "point3") {
      value.type = pbrt_type::point;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "normal" || type == "normal3") {
      value.type = pbrt_type::normal;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "vector" || type == "vector3") {
      value.type = pbrt_type::vector;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "point2") {
      value.type = pbrt_type::point2;
      if (!parse_pvalues(str, value.value2f, value.vector2f)) return false;
    } else if (type == "vector2") {
      value.type = pbrt_type::vector2;
      if (!parse_pvalues(str, value.value2f, value.vector2f)) return false;
    } else if (type == "blackbody") {
      value.type = pbrt_type::color;
      // auto blackbody = vec2f{0, 0};
      // auto vec tor2f  = vector<array<float, 2>>{};
      // parse_pvalues(str, blackbody, vector2f);
      // if (!vector2f.empty()) return false;
      // value.value3f = blackbody_to_rgb(blackbody[0]) * blackbody[1];
      auto blackbody = 0.0f;
      auto vector1f  = vector<float>{};
      if (!parse_pvalues(str, blackbody, vector1f)) return false;
      if (vector1f.size() < 2) {
        value.value3f = blackbody_to_rgb_(blackbody);
      } else {
        value.value3f = mul(blackbody_to_rgb_(vector1f[0]), vector1f[1]);
      }
    } else if (type == "color" || type == "rgb") {
      value.type = pbrt_type::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "xyz") {
      value.type = pbrt_type::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
      // xyz conversion
      return false;
    } else if (type == "spectrum") {
      auto is_string = false;
      auto str1      = str;
      skip_whitespace(str1);
      if (!str1.empty() && str1.front() == '"') {
        is_string = true;
      } else if (!str1.empty() && str1.front() == '[') {
        str1.remove_prefix(1);
        skip_whitespace(str1);
        if (!str1.empty() && str1.front() == '"') is_string = true;
      }
      if (is_string) {
        value.type     = pbrt_type::color;
        auto filename  = ""s;
        auto filenames = vector<string>{};
        skip_whitespace(str);
        auto has_parens = str.front() == '[';
        if (has_parens) str.remove_prefix(1);
        if (!parse_value(str, filename)) return false;
        if (has_parens) {
          skip_whitespace(str);
          if (str.front() != ']') return false;
          str.remove_prefix(1);
        }
        if (str.empty()) return false;
        auto filenamep = path_filename(filename);
        auto name      = string_view{filenamep};
        if (ends_with(name, ".spd")) {
          name.remove_suffix(4);
          if (name == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (ends_with(name, ".eta")) {
            name.remove_suffix(4);
            auto eta      = get_etak(string{name}).first;
            value.value3f = {eta[0], eta[1], eta[2]};
          } else if (ends_with(name, ".k")) {
            name.remove_suffix(2);
            auto k        = get_etak(string{name}).second;
            value.value3f = {k[0], k[1], k[2]};
          } else {
            return false;
          }
        } else if (starts_with(name, "metal-")) {
          name.remove_prefix(6);
          if (ends_with(name, "-eta")) {
            name.remove_suffix(4);
            auto eta      = get_etak(string{name}).first;
            value.value3f = {eta[0], eta[1], eta[2]};
          } else if (ends_with(name, "-k")) {
            name.remove_suffix(2);
            auto k        = get_etak(string{name}).second;
            value.value3f = {k[0], k[1], k[2]};
          } else {
            return false;
          }
        } else if (starts_with(name, "glass-")) {
          value.value3f = {1.5, 1.5, 1.5};
        } else {
          return false;
        }
      } else {
        value.type = pbrt_type::spectrum;
        if (!parse_pvalues(str, value.value1f, value.vector1f)) return false;
      }
    } else {
      return false;
    }
    skip_whitespace(str);
  }
  return true;
}

// Other pbrt elements
struct pbrt_film {
  // film approximation
  string        filename   = "";
  array<int, 2> resolution = {0, 0};
};

// Pbrt area light
struct pbrt_arealight {
  // arealight parameters
  string          name     = "";
  array<float, 3> emission = {0, 0, 0};
};

// Pbrt medium. Not parsed at the moment.
struct pbrt_medium {
  // medium parameters
  string name = "";
};

// convert pbrt films
[[nodiscard]] static bool convert_film(pbrt_film& film,
    const pbrt_command& command, const string& filename, bool verbose = false) {
  if (command.type == "image") {
    film.resolution = {512, 512};
    get_pbrt_value(command.values, "xresolution", film.resolution[0]);
    get_pbrt_value(command.values, "yresolution", film.resolution[1]);
    film.filename = "out.png"s;
    get_pbrt_value(command.values, "filename", film.filename);
  } else if (command.type == "rgb") {
    film.resolution = {512, 512};
    get_pbrt_value(command.values, "xresolution", film.resolution[0]);
    get_pbrt_value(command.values, "yresolution", film.resolution[1]);
    film.filename = "out.png"s;
    get_pbrt_value(command.values, "filename", film.filename);
  } else {
    return false;
  }
  return true;
}

// convert pbrt elements
[[nodiscard]] static bool convert_camera(pbrt_camera& pcamera,
    const pbrt_command& command, const array<int, 2>& resolution,
    const string& filename, bool verbose = false) {
  auto cframe        = inverse(command.frame);
  cframe[2]          = neg(cframe[2]);
  auto cfrend        = inverse(command.frend);
  cfrend[2]          = neg(cfrend[2]);
  pcamera.frame      = flatten(cframe);
  pcamera.frend      = flatten(cfrend);
  pcamera.resolution = resolution;
  auto film_aspect   = (resolution[0] == 0 || resolution[1] == 0)
                           ? 1
                           : (float)resolution[0] / (float)resolution[1];
  if (command.type == "perspective") {
    auto fov = 90.0f;
    get_pbrt_value(command.values, "fov", fov);
    // auto lensradius = if(!get_pbrt_value(values, "lensradius", 0.0f);
    pcamera.aspect = film_aspect;
    if (pcamera.aspect >= 1) {
      pcamera.lens = (0.036f / pcamera.aspect) /
                     (2 * std::tan((fov * (float)M_PI / 180) / 2));
    } else {
      pcamera.lens = (0.036f * pcamera.aspect) /
                     (2 * std::tan((fov * (float)M_PI / 180) / 2));
    }
    get_pbrt_value(command.values, "frameaspectratio", pcamera.aspect);
    pcamera.focus = 10.0f;
    get_pbrt_value(command.values, "focaldistance", pcamera.focus);
  } else if (command.type == "realistic") {
    auto lensfile = ""s;
    get_pbrt_value(command.values, "lensfile", lensfile);
    lensfile     = lensfile.substr(0, lensfile.size() - 4);
    lensfile     = lensfile.substr(lensfile.find('.') + 1);
    lensfile     = lensfile.substr(0, lensfile.size() - 2);
    auto lens    = std::max((float)std::atof(lensfile.c_str()), 35.0f) * 0.001f;
    pcamera.lens = 2 * atan(0.036f / (2 * lens));
    pcamera.aperture = 0.0f;
    get_pbrt_value(command.values, "aperturediameter", pcamera.aperture);
    pcamera.focus = 10.0f;
    get_pbrt_value(command.values, "focusdistance", pcamera.focus);
    pcamera.aspect = film_aspect;
  } else {
    return false;
  }
  return true;
}

// convert pbrt textures
[[nodiscard]] static bool convert_texture(pbrt_texture& ptexture,
    const pbrt_command&                                 command,
    unordered_map<string, pbrt_texture>& texture_map, const string& filename,
    bool verbose = false) {
  auto make_filename = [&texture_map](const string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return pos->second.filename;
  };

  ptexture.name = command.name;
  if (command.type == "imagemap") {
    ptexture.filename = "";
    get_pbrt_value(command.values, "filename", ptexture.filename);
  } else if (command.type == "constant") {
    ptexture.constant = array<float, 3>{1, 1, 1};
    get_pbrt_value(command.values, "value", ptexture.constant);
  } else if (command.type == "bilerp") {
    ptexture.constant = {1, 0, 0};
  } else if (command.type == "checkerboard") {
    // auto tex1     = if(!get_pbrt_value(command.values, "tex1",
    // pair{array<float, 3>{1,1,1},
    // ""s}); auto tex2     = if(!get_pbrt_value(command.values, "tex2",
    //  pair{array<float, 3>{0}, ""s}); auto rgb1     = tex1.second == "" ?
    //  tex1.first :
    // array<float, 3>{0.4f, 0.4f, 0.4f}; auto rgb2     = tex1.second == "" ?
    // tex2.first : array<float, 3>{0.6f, 0.6f, 0.6f}; auto params   =
    // proc_image_params{}; params.type = proc_image_params::type_t::checker;
    // params.color0 = {rgb1[0], rgb1[1], rgb1[2], 1}; params.color1 = {rgb2[0],
    // rgb2[1], rgb2[2], 1}; params.scale = 2; make_proc_image(texture.hdr,
    // params); float_to_byte(texture.ldr, texture.hdr); texture.hdr = {};
    ptexture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "dots") {
    ptexture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "fbm") {
    ptexture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "marble") {
    ptexture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "mix") {
    auto tex1 = pair{array<float, 3>{0, 0, 0}, ""s},
         tex2 = pair{array<float, 3>{1, 1, 1}, ""s};
    get_pbrt_value(command.values, "tex1", tex1);
    get_pbrt_value(command.values, "tex2", tex2);
    if (!make_filename(tex1.second).empty()) {
      ptexture.filename = make_filename(tex1.second);
    } else if (!make_filename(tex2.second).empty()) {
      ptexture.filename = make_filename(tex2.second);
    } else {
      ptexture.constant = {1, 0, 0};
    }
  } else if (command.type == "scale") {
    auto tex1 = pair{array<float, 3>{1, 1, 1}, ""s},
         tex2 = pair{array<float, 3>{1, 1, 1}, ""s};
    get_pbrt_value(command.values, "tex1", tex2);
    get_pbrt_value(command.values, "tex2", tex1);
    if (!make_filename(tex1.second).empty()) {
      ptexture.filename = make_filename(tex1.second);
    } else if (!make_filename(tex2.second).empty()) {
      ptexture.filename = make_filename(tex2.second);
    } else {
      ptexture.constant = {1, 0, 0};
    }
  } else if (command.type == "uv") {
    ptexture.constant = {1, 0, 0};
  } else if (command.type == "windy") {
    ptexture.constant = {1, 0, 0};
  } else if (command.type == "wrinkled") {
    ptexture.constant = {1, 0, 0};
  } else {
    return false;
  }
  return true;
}

// convert pbrt materials
[[nodiscard]] static bool convert_material(pbrt_material& pmaterial,
    const pbrt_command& command, unordered_map<string, int>& texture_map,
    const unordered_map<string, pbrt_material>& named_materials,
    const unordered_map<string, pbrt_texture>&  named_textures,
    const string& filename, bool verbose = false) {
  // helpers
  auto get_texture_id = [&texture_map](const string& path) {
    if (path.empty()) return -1;
    auto texture_it = texture_map.find(path);
    if (texture_it == texture_map.end()) {
      auto texture_id   = (int)texture_map.size();
      texture_map[path] = texture_id;
      return texture_id;
    } else {
      return texture_it->second;
    }
  };
  auto get_texture = [&](const vector<pbrt_value>& values, const string& name,
                         array<float, 3>& color, int& texture_id,
                         const array<float, 3>& def) {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      color      = textured.first;
      texture_id = -1;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color      = texture.constant;
        texture_id = -1;
      } else {
        color      = {1, 1, 1};
        texture_id = get_texture_id(texture.filename);
      }
    }
  };
  auto get_scalar = [&](const vector<pbrt_value>& values, const string& name,
                        float& scalar, float def) {
    auto textured = pair{array<float, 3>{def, def, def}, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      scalar = (textured.first[0] + textured.first[1] + textured.first[2]) / 3;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        scalar =
            (texture.constant[0] + texture.constant[1] + texture.constant[2]) /
            3;
      } else {
        scalar = def;
      }
    }
  };
  auto get_color = [&](const vector<pbrt_value>& values, const string& name,
                       array<float, 3>& color, const array<float, 3>& def) {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      color = textured.first;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color = texture.constant;
      } else {
        color = def;
      }
    }
  };

  auto get_roughness = [&](const vector<pbrt_value>& values, float& roughness,
                           float def = 0.1f) {
    auto roughness_ = pair{array<float, 3>{def, def, def}, ""s};
    get_pbrt_value(values, "roughness", roughness_);
    auto uroughness = roughness_, vroughness = roughness_;
    auto remaproughness = true;
    get_pbrt_value(values, "uroughness", uroughness);
    get_pbrt_value(values, "vroughness", vroughness);
    get_pbrt_value(values, "remaproughness", remaproughness);

    roughness = 0;
    if (uroughness.first == array<float, 3>{0, 0, 0} ||
        vroughness.first == array<float, 3>{0, 0, 0})
      return;
    auto uroughness_mean =
        (uroughness.first[0] + uroughness.first[1] + uroughness.first[2]) / 3;
    auto vroughness_mean =
        (vroughness.first[0] + vroughness.first[1] + vroughness.first[2]) / 3;
    roughness = (uroughness_mean + vroughness_mean) / 2;
    // from pbrt code
    if (remaproughness) {
      roughness = std::max(roughness, 1e-3f);
      auto x    = log(roughness);
      roughness = 1.62142f + 0.819955f * x + 0.1734f * x * x +
                  0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }
    roughness = sqrt(roughness);
  };

  auto eta_to_reflectivity1 = [](float eta, float etak) -> float {
    return ((eta - 1) * (eta - 1) + etak * etak) /
           ((eta + 1) * (eta + 1) + etak * etak);
  };
  auto eta_to_reflectivity =
      [&](const array<float, 3>& eta,
          const array<float, 3>& etak = {0, 0, 0}) -> array<float, 3> {
    return {eta_to_reflectivity1(eta[0], etak[0]),
        eta_to_reflectivity1(eta[1], etak[1]),
        eta_to_reflectivity1(eta[2], etak[2])};
  };

  pmaterial.name = command.name;
  if (command.type == "uber") {
    auto diffuse      = array<float, 3>{0, 0, 0},
         specular     = array<float, 3>{0, 0, 0},
         transmission = array<float, 3>{0, 0, 0};
    auto diffuse_map = -1, specular_map = -1, transmission_map = -1;
    get_texture(command.values, "Kd", diffuse, diffuse_map,
        array<float, 3>{0.25, 0.25, 0.25});
    get_texture(command.values, "Ks", specular, specular_map,
        array<float, 3>{0.25, 0.25, 0.25});
    get_texture(command.values, "Kt", transmission, transmission_map,
        array<float, 3>{0, 0, 0});
    if (std::max(transmission[0], std::max(transmission[1], transmission[2])) >
        0.1) {
      pmaterial.type      = pbrt_mtype::thinglass;
      pmaterial.color     = transmission;
      pmaterial.color_tex = transmission_map;
    } else if (std::max(specular[0], std::max(specular[1], specular[2])) >
               0.1) {
      pmaterial.type      = pbrt_mtype::plastic;
      pmaterial.color     = diffuse;
      pmaterial.color_tex = diffuse_map;
    } else {
      pmaterial.type      = pbrt_mtype::plastic;
      pmaterial.color     = diffuse;
      pmaterial.color_tex = diffuse_map;
    }
    get_scalar(command.values, "opacity", pmaterial.opacity, 1);
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    get_roughness(command.values, pmaterial.roughness, 0.1f);
  } else if (command.type == "plastic") {
    pmaterial.type = pbrt_mtype::plastic;
    get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.25, 0.25, 0.25});
    // get_scalar(command.values, "Ks", pmaterial.specular, 0.25))
    //   return parse_error();
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0.1f;
    get_roughness(command.values, pmaterial.roughness, 0.1f);
  } else if (command.type == "coateddiffuse") {
    pmaterial.type = pbrt_mtype::plastic;
    get_texture(command.values, "reflectance", pmaterial.color,
        pmaterial.color_tex, array<float, 3>{0.25, 0.25, 0.25});
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0.1f;
    get_roughness(command.values, pmaterial.roughness, 0.1f);
  } else if (command.type == "translucent") {
    // not well supported yet
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.25, 0.25, 0.25});
    // get_scalar(command.values, "Ks", pmaterial.specular, 0.25))
    //   return parse_error();
    // get_scalar(command.values, "eta", pmaterial.ior, 1.5))
    //   return parse_error();
    // get_roughness(command.values, pmaterial.roughness, 0.1))
    //   return parse_error();
  } else if (command.type == "diffusetransmission") {
    // not well supported yet
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "reflectance", pmaterial.color,
        pmaterial.color_tex, array<float, 3>{0.25f, 0.25f, 0.25f});
    // get_texture(command.values, "transmittance", pmaterial.color,
    //         pmaterial.color_tex, array<float, 3>{0.25, 0.25, 0.25}))
    //   return parse_error();
  } else if (command.type == "matte") {
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.5, 0.5, 0.5});
  } else if (command.type == "diffuse") {
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "reflectance", pmaterial.color,
        pmaterial.color_tex, array<float, 3>{0.5f, 0.5f, 0.5f});
  } else if (command.type == "mirror") {
    pmaterial.type = pbrt_mtype::metal;
    get_texture(command.values, "Kr", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.9f, 0.9f, 0.9f});
    pmaterial.roughness = 0;
  } else if (command.type == "metal") {
    pmaterial.type = pbrt_mtype::metal;
    // get_texture(
    //     values, "Kr", material->specular, material->specular_tex,
    //     array<float, 3>{1,1,1});
    auto eta = array<float, 3>{0, 0, 0}, etak = array<float, 3>{0, 0, 0};
    get_color(command.values, "eta", eta,
        array<float, 3>{0.2004376970f, 0.9240334304f, 1.1022119527f});
    get_color(command.values, "k", etak,
        array<float, 3>{3.9129485033f, 2.4528477015f, 2.1421879552f});
    pmaterial.color     = eta_to_reflectivity(eta, etak);
    pmaterial.roughness = 0.01f;
    get_roughness(command.values, pmaterial.roughness, 0.01f);
  } else if (command.type == "conductor") {
    pmaterial.type = pbrt_mtype::metal;
    auto eta = array<float, 3>{0, 0, 0}, etak = array<float, 3>{0, 0, 0};
    get_color(command.values, "eta", eta,
        array<float, 3>{0.2004376970f, 0.9240334304f, 1.1022119527f});
    get_color(command.values, "k", etak,
        array<float, 3>{3.9129485033f, 2.4528477015f, 2.1421879552f});
    pmaterial.color     = eta_to_reflectivity(eta, etak);
    pmaterial.roughness = 0.01f;
    get_roughness(command.values, pmaterial.roughness, 0.01f);
  } else if (command.type == "coatedconductor") {
    pmaterial.type = pbrt_mtype::metal;
    auto eta = array<float, 3>{0, 0, 0}, etak = array<float, 3>{0, 0, 0};
    get_color(command.values, "conductor.eta", eta,
        array<float, 3>{0.2004376970f, 0.9240334304f, 1.1022119527f});
    get_color(command.values, "conductor.k", etak,
        array<float, 3>{3.9129485033f, 2.4528477015f, 2.1421879552f});
    pmaterial.color     = eta_to_reflectivity(eta, etak);
    pmaterial.roughness = 0.01f;
    get_roughness(command.values, pmaterial.roughness, 0.01f);
  } else if (command.type == "substrate") {
    // not well supported
    pmaterial.type = pbrt_mtype::plastic;
    get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.5f, 0.5f, 0.5f});
    auto specular = 0.0f;
    get_scalar(command.values, "Ks", specular, 0.5f);
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0.1f;
    get_roughness(command.values, pmaterial.roughness, 0.1f);
  } else if (command.type == "glass") {
    pmaterial.type = pbrt_mtype::glass;
    get_texture(command.values, "Kt", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{1, 1, 1});
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0;
    get_roughness(command.values, pmaterial.roughness, 0.0f);
  } else if (command.type == "dielectric") {
    pmaterial.type  = pbrt_mtype::glass;
    pmaterial.color = {1, 1, 1};
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0;
    get_roughness(command.values, pmaterial.roughness, 0.0f);
  } else if (command.type == "thindielectric") {
    pmaterial.type  = pbrt_mtype::thinglass;
    pmaterial.color = {1, 1, 1};
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0;
    get_roughness(command.values, pmaterial.roughness, 0.0f);
  } else if (command.type == "hair") {
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "color", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0, 0, 0});
    pmaterial.roughness = 1;
    if (verbose) printf("hair material not properly supported\n");
  } else if (command.type == "disney") {
    pmaterial.type = pbrt_mtype::matte;
    get_texture(command.values, "color", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.5f, 0.5f, 0.5f});
    pmaterial.roughness = 1;
    if (verbose) printf("disney material not properly supported\n");
  } else if (command.type == "kdsubsurface") {
    pmaterial.type = pbrt_mtype::plastic;
    get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
        array<float, 3>{0.5f, 0.5f, 0.5f});
    // get_scalar(command.values, "Kr", pmaterial.specular, 1))
    //   return parse_error();
    get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
    pmaterial.roughness = 0;
    get_roughness(command.values, pmaterial.roughness, 0);
    if (verbose) printf("kdsubsurface material not properly supported\n");
  } else if (command.type == "subsurface") {
    pmaterial.type = pbrt_mtype::subsurface;
    // get_scalar(command.values, "Kr", pmaterial.specular, 1))
    //   return parse_error();
    // get_scalar(command.values, "Kt", pmaterial.transmission, 1))
    //   return parse_error();
    pmaterial.color = {1, 1, 1};
    get_scalar(command.values, "eta", pmaterial.ior, 1.5);
    pmaterial.roughness = 0;
    get_roughness(command.values, pmaterial.roughness, 0);
    auto scale = 1.0f;
    get_pbrt_value(command.values, "scale", scale);
    pmaterial.volscale = 1 / scale;
    auto sigma_a = array<float, 3>{0, 0, 0}, sigma_s = array<float, 3>{0, 0, 0};
    auto sigma_a_tex = -1, sigma_s_tex = -1;
    get_texture(command.values, "sigma_a", sigma_a, sigma_a_tex,
        array<float, 3>{0.011f, .0024f, .014f});
    get_texture(command.values, "sigma_prime_s", sigma_s, sigma_s_tex,
        array<float, 3>{2.55f, 3.12f, 3.77f});
    pmaterial.volmeanfreepath = {1 / (sigma_a[0] + sigma_s[0]),
        1 / (sigma_a[1] + sigma_s[1]), 1 / (sigma_a[2] + sigma_s[2])};
    pmaterial.volscatter      = {sigma_s[0] / (sigma_a[0] + sigma_s[0]),
        sigma_s[1] / (sigma_a[1] + sigma_s[1]),
        sigma_s[2] / (sigma_a[2] + sigma_s[2])};
    if (verbose) printf("subsurface material not properly supported\n");
  } else if (command.type == "mix") {
    auto namedmaterial1 = ""s, namedmaterial2 = ""s;
    get_pbrt_value(command.values, "namedmaterial1", namedmaterial1);
    get_pbrt_value(command.values, "namedmaterial2", namedmaterial2);
    auto matname = (!namedmaterial1.empty()) ? namedmaterial1 : namedmaterial2;
    auto matit   = named_materials.find(matname);
    if (matit == named_materials.end()) return false;
    auto saved_name = pmaterial.name;
    pmaterial       = matit->second;
    pmaterial.name  = saved_name;
    if (verbose) printf("mix material not properly supported\n");
  } else if (command.type == "fourier") {
    auto bsdffile = ""s;
    get_pbrt_value(command.values, "bsdffile", bsdffile);
    if (bsdffile.rfind('/') != string::npos)
      bsdffile = bsdffile.substr(bsdffile.rfind('/') + 1);
    if (bsdffile == "paint.bsdf") {
      pmaterial.type      = pbrt_mtype::plastic;
      pmaterial.color     = {0.6f, 0.6f, 0.6f};
      pmaterial.ior       = 1.5f;
      pmaterial.roughness = 0.2f;
    } else if (bsdffile == "ceramic.bsdf") {
      pmaterial.type      = pbrt_mtype::plastic;
      pmaterial.color     = {0.6f, 0.6f, 0.6f};
      pmaterial.ior       = 1.5f;
      pmaterial.roughness = 0.25f;
    } else if (bsdffile == "leather.bsdf") {
      pmaterial.type      = pbrt_mtype::plastic;
      pmaterial.color     = {0.6f, 0.57f, 0.48f};
      pmaterial.ior       = 1.5f;
      pmaterial.roughness = 0.3f;
    } else if (bsdffile == "coated_copper.bsdf") {
      pmaterial.type = pbrt_mtype::metal;
      auto eta  = array<float, 3>{0.2004376970f, 0.9240334304f, 1.1022119527f};
      auto etak = array<float, 3>{3.9129485033f, 2.4528477015f, 2.1421879552f};
      pmaterial.color     = eta_to_reflectivity(eta, etak);
      pmaterial.roughness = 0.01f;
    } else if (bsdffile == "roughglass_alpha_0.2.bsdf") {
      pmaterial.type      = pbrt_mtype::glass;
      pmaterial.color     = {1, 1, 1};
      pmaterial.ior       = 1.5f;
      pmaterial.roughness = 0.2f;
    } else if (bsdffile == "roughgold_alpha_0.2.bsdf") {
      pmaterial.type = pbrt_mtype::metal;
      auto eta  = array<float, 3>{0.1431189557f, 0.3749570432f, 1.4424785571f};
      auto etak = array<float, 3>{3.9831604247f, 2.3857207478f, 1.6032152899f};
      pmaterial.color     = eta_to_reflectivity(eta, etak);
      pmaterial.roughness = 0.2f;
    } else {
      return false;
    }
  } else {
    return false;
  }
  return true;
}

// Make a triangle shape from a quad grid
template <typename PositionFunc, typename NormalFunc>
static void make_shape(vector<array<int, 3>>& triangles,
    vector<array<float, 3>>& positions, vector<array<float, 3>>& normals,
    vector<array<float, 2>>& texcoords, const array<int, 2>& steps,
    const PositionFunc& position_func, const NormalFunc& normal_func) {
  auto vid = [steps](int i, int j) { return j * (steps[0] + 1) + i; };
  auto tid = [steps](
                 int i, int j, int c) { return (j * steps[0] + i) * 2 + c; };
  positions.resize((steps[0] + 1) * (steps[1] + 1));
  normals.resize((steps[0] + 1) * (steps[1] + 1));
  texcoords.resize((steps[0] + 1) * (steps[1] + 1));
  for (auto j = 0; j < steps[1] + 1; j++) {
    for (auto i = 0; i < steps[0] + 1; i++) {
      auto uv = array<float, 2>{i / (float)steps[0], j / (float)steps[1]};
      positions[vid(i, j)] = position_func(uv);
      normals[vid(i, j)]   = normal_func(uv);
      texcoords[vid(i, j)] = uv;
    }
  }
  triangles.resize(steps[0] * steps[1] * 2);
  for (auto j = 0; j < steps[1]; j++) {
    for (auto i = 0; i < steps[0]; i++) {
      triangles[tid(i, j, 0)] = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
      triangles[tid(i, j, 1)] = {vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)};
    }
  }
}

// pbrt sphere
static void make_sphere(vector<array<int, 3>>& triangles,
    vector<array<float, 3>>& positions, vector<array<float, 3>>& normals,
    vector<array<float, 2>>& texcoords, const array<int, 2>& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const array<float, 2>& uv) {
        auto pt = array<float, 2>{
            2 * (float)M_PI * uv[0], (float)M_PI * (1 - uv[1])};
        return array<float, 3>{radius * cos(pt[0]) * sin(pt[1]),
            radius * sin(pt[0]) * sin(pt[1]), radius * cos(pt[1])};
      },
      [](const array<float, 2>& uv) {
        auto pt = array<float, 2>{
            2 * (float)M_PI * uv[0], (float)M_PI * (1 - uv[1])};
        return array<float, 3>{
            cos(pt[0]) * sin(pt[1]), sin(pt[0]) * sin(pt[1]), cos(pt[1])};
      });
}
static void make_disk(vector<array<int, 3>>& triangles,
    vector<array<float, 3>>& positions, vector<array<float, 3>>& normals,
    vector<array<float, 2>>& texcoords, const array<int, 2>& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const array<float, 2>& uv) {
        auto a = 2 * (float)M_PI * uv[0];
        return array<float, 3>{
            radius * (1 - uv[1]) * cos(a), radius * (1 - uv[1]) * sin(a), 0};
      },
      [](const array<float, 2>& uv) {
        return array<float, 3>{0, 0, 1};
      });
}
static void make_quad(vector<array<int, 3>>& triangles,
    vector<array<float, 3>>& positions, vector<array<float, 3>>& normals,
    vector<array<float, 2>>& texcoords, const array<int, 2>& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const array<float, 2>& uv) {
        return array<float, 3>{
            (uv[0] - 0.5f) * radius, (uv[1] - 0.5f) * radius, 0};
      },
      [](const array<float, 2>& uv) {
        return array<float, 3>{0, 0, 1};
      });
}

// Convert pbrt shapes
[[nodiscard]] static bool convert_shape(pbrt_shape& pshape,
    const pbrt_command& command, string& alphamap,
    const unordered_map<string, pbrt_texture>& named_textures,
    const string& ply_dirname, bool ply_meshes, const string& filename,
    bool verbose = false) {
  // helpers
  auto get_alpha = [&](const vector<pbrt_value>& values, const string& name,
                       string& filename) -> bool {
    auto def      = 1.0f;
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      filename = "";
    } else {
      filename = named_textures.at(textured.second).filename;
    }
    return true;
  };

  pshape.frame = flatten(command.frame);
  pshape.frend = flatten(command.frend);
  if (command.type == "trianglemesh") {
    pshape.positions = {};
    pshape.normals   = {};
    pshape.texcoords = {};
    pshape.triangles = {};
    get_pbrt_value(command.values, "P", pshape.positions);
    get_pbrt_value(command.values, "N", pshape.normals);
    get_pbrt_value(command.values, "uv", pshape.texcoords);
    for (auto& uv : pshape.texcoords) uv[1] = (1 - uv[1]);
    get_pbrt_value(command.values, "indices", pshape.triangles);
  } else if (command.type == "loopsubdiv") {
    pshape.positions = {};
    pshape.triangles = {};
    get_pbrt_value(command.values, "P", pshape.positions);
    get_pbrt_value(command.values, "indices", pshape.triangles);
    pshape.normals.resize(pshape.positions.size());
    // compute_normals(pshape.normals, pshape.triangles, pshape.positions);
  } else if (command.type == "plymesh") {
    pshape.filename_ = ""s;
    get_pbrt_value(command.values, "filename", pshape.filename_);
    get_alpha(command.values, "alpha", alphamap);
    if (ply_meshes) {
      auto error = string{};
      auto ply   = ply_model{};
      if (!load_ply(path_join(ply_dirname, pshape.filename_), ply, error))
        return false;
      get_positions(ply, pshape.positions);
      get_normals(ply, pshape.normals);
      get_texcoords(ply, pshape.texcoords);
      get_triangles(ply, pshape.triangles);
    }
  } else if (command.type == "sphere") {
    auto radius = 1.0f;
    get_pbrt_value(command.values, "radius", radius);
    make_sphere(pshape.triangles, pshape.positions, pshape.normals,
        pshape.texcoords, {32, 16}, radius);
  } else if (command.type == "disk") {
    auto radius = 1.0f;
    get_pbrt_value(command.values, "radius", radius);
    make_disk(pshape.triangles, pshape.positions, pshape.normals,
        pshape.texcoords, {32, 1}, radius);
  } else {
    return false;
  }
  return true;
}

// Convert pbrt arealights
[[nodiscard]] static bool convert_arealight(pbrt_arealight& parealight,
    const pbrt_command& command, const string& filename, bool verbose = false) {
  parealight.name = command.name;
  if (command.type == "diffuse") {
    auto l = array<float, 3>{1, 1, 1}, scale = array<float, 3>{1, 1, 1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    parealight.emission = {l[0] * scale[0], l[1] * scale[1], l[2] * scale[2]};
  } else {
    return false;
  }
  return true;
}

// Convert pbrt lights
[[nodiscard]] static bool convert_light(pbrt_light& plight,
    const pbrt_command& command, const string& filename, bool verbose = false) {
  plight.frame = flatten(command.frame);
  plight.frend = flatten(command.frend);
  if (command.type == "distant") {
    auto l = array<float, 3>{1, 1, 1}, scale = array<float, 3>{1, 1, 1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    plight.emission = {l[0] * scale[0], l[1] * scale[1], l[2] * scale[2]};
    plight.from     = array<float, 3>{0, 0, 0};
    plight.to       = array<float, 3>{0, 0, 1};
    get_pbrt_value(command.values, "from", plight.from);
    get_pbrt_value(command.values, "to", plight.to);
    plight.distant       = true;
    auto distant_dist    = 100.0f;
    auto size            = distant_dist * std::sin(5 * (float)M_PI / 180);
    auto dscale          = (distant_dist * distant_dist) / (size * size);
    plight.area_emission = {plight.emission[0] * dscale,
        plight.emission[1] * dscale, plight.emission[2] * dscale};
    plight.area_frame    = flatten(mul(unflatten(plight.frame),
           lookat_frame(mul(normalize(sub(plight.from, plight.to)), distant_dist),
               {0, 0, 0}, {0, 1, 0}, true)));
    plight.area_frend    = flatten(mul(unflatten(plight.frend),
           lookat_frame(mul(normalize(sub(plight.from, plight.to)), distant_dist),
               {0, 0, 0}, {0, 1, 0}, true)));
    auto texcoords       = vector<array<float, 2>>{};
    make_quad(plight.area_triangles, plight.area_positions, plight.area_normals,
        texcoords, {4, 2}, size);
  } else if (command.type == "point" || command.type == "goniometric" ||
             command.type == "spot") {
    auto i = array<float, 3>{1, 1, 1}, scale = array<float, 3>{1, 1, 1};
    get_pbrt_value(command.values, "I", i);
    get_pbrt_value(command.values, "scale", scale);
    plight.emission = {i[0] * scale[0], i[1] * scale[1], i[2] * scale[2]};
    plight.from     = {0, 0, 0};
    get_pbrt_value(command.values, "from", plight.from);
    plight.area_emission = plight.emission;
    plight.area_frame    = flatten(
           mul(unflatten(plight.frame), translation_frame(plight.from)));
    plight.area_frend = flatten(
        mul(unflatten(plight.frend), translation_frame(plight.from)));
    auto texcoords = vector<array<float, 2>>{};
    make_sphere(plight.area_triangles, plight.area_positions,
        plight.area_normals, texcoords, {4, 2}, 0.0025f);
  } else {
    return false;
  }
  return true;
}

[[nodiscard]] static bool convert_environment(pbrt_environment& penvironment,
    const pbrt_command& command, unordered_map<string, int>& texture_map,
    const string& filename, bool verbose = false) {
  penvironment.frame = flatten(command.frame);
  penvironment.frend = flatten(command.frend);
  penvironment.frame = flatten(mul(unflatten(penvironment.frame),
      array<array<float, 3>, 4>{array<float, 3>{1, 0, 0},
          array<float, 3>{0, 0, 1}, array<float, 3>{0, 1, 0},
          array<float, 3>{0, 0, 0}}));
  penvironment.frend = flatten(mul(unflatten(penvironment.frend),
      array<array<float, 3>, 4>{array<float, 3>{1, 0, 0},
          array<float, 3>{0, 0, 1}, array<float, 3>{0, 1, 0},
          array<float, 3>{0, 0, 0}}));
  if (command.type == "infinite") {
    auto l = array<float, 3>{1, 1, 1}, scale = array<float, 3>{1, 1, 1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    penvironment.emission = {l[0] * scale[0], l[1] * scale[1], l[2] * scale[2]};
    penvironment.emission_tex = -1;
    auto mapname              = ""s;
    get_pbrt_value(command.values, "mapname", mapname);
    if (!mapname.empty()) {
      if (texture_map.find(mapname) == texture_map.end()) {
        auto texture_id      = (int)texture_map.size();
        texture_map[mapname] = texture_id;
      }
      penvironment.emission_tex = texture_map.at(mapname);
    }
  } else {
    return false;
  }
  return true;
}

// pbrt stack ctm
struct pbrt_stack_element {
  array<array<float, 3>, 4> transform_start        = {array<float, 3>{1, 0, 0},
      array<float, 3>{0, 1, 0}, array<float, 3>{0, 0, 1},
      array<float, 3>{0, 0, 0}};
  array<array<float, 3>, 4> transform_end          = {array<float, 3>{1, 0, 0},
      array<float, 3>{0, 1, 0}, array<float, 3>{0, 0, 1},
      array<float, 3>{0, 0, 0}};
  pbrt_material             material               = {};
  pbrt_arealight            arealight              = {};
  pbrt_medium               interior               = {};
  pbrt_medium               exterior               = {};
  bool                      reverse                = false;
  bool                      active_transform_start = true;
  bool                      active_transform_end   = true;
};

// pbrt parsing context
struct pbrt_context {
  vector<pbrt_stack_element>                stack           = {};
  unordered_map<string, pbrt_stack_element> coordsys        = {};
  string                                    cur_object      = "";
  array<int, 2>                             film_resolution = {512, 512};
};

// load pbrt
static bool load_pbrt(const string& filename, pbrt_model& pbrt, string& error,
    pbrt_context& ctx, unordered_map<string, int>& material_map,
    unordered_map<string, int>&           texture_map,
    unordered_map<string, pbrt_material>& named_materials,
    unordered_map<string, pbrt_texture>&  named_textures,
    unordered_map<string, pbrt_medium>&   named_mediums,
    unordered_map<string, vector<int>>&   named_objects,
    const string& ply_dirname, bool ply_meshes) {
  // load data
  auto data = string{};
  if (!load_text(filename, data, error)) return false;

  // helpers
  auto set_transform = [](pbrt_stack_element&               ctx,
                           const array<array<float, 3>, 4>& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](pbrt_stack_element&               ctx,
                              const array<array<float, 3>, 4>& xform) {
    if (ctx.active_transform_start)
      ctx.transform_start = mul(ctx.transform_start, xform);
    if (ctx.active_transform_end)
      ctx.transform_end = mul(ctx.transform_end, xform);
  };

  // init stack
  if (ctx.stack.empty()) ctx.stack.emplace_back();

  // parse command by command
  auto data_view   = string_view{data.data(), data.size()};
  auto line        = ""s;
  auto parse_error = [&filename, &error]() {
    error = "cannot parse " + filename;
    return false;
  };
  auto dependent_error = [&filename, &error]() {
    error = "cannot load " + filename + " since " + error;
    return false;
  };
  while (read_pbrt_cmdline(data_view, line)) {
    auto str = string_view{line};
    // get command
    auto cmd = ""s;
    if (!parse_command(str, cmd)) return parse_error();
    if (cmd == "WorldBegin") {
      ctx.stack.push_back({});
    } else if (cmd == "WorldEnd") {
      if (ctx.stack.empty()) return parse_error();
      ctx.stack.pop_back();
      if (ctx.stack.size() != 1) return parse_error();
    } else if (cmd == "AttributeBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "AttributeEnd") {
      if (ctx.stack.empty()) return parse_error();
      ctx.stack.pop_back();
    } else if (cmd == "TransformBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "TransformEnd") {
      if (ctx.stack.empty()) return parse_error();
      ctx.stack.pop_back();
    } else if (cmd == "ObjectBegin") {
      ctx.stack.push_back(ctx.stack.back());
      if (!parse_param(str, ctx.cur_object)) return parse_error();
      named_objects[ctx.cur_object] = {};
    } else if (cmd == "ObjectEnd") {
      ctx.stack.pop_back();
      ctx.cur_object = "";
    } else if (cmd == "ObjectInstance") {
      auto object = ""s;
      if (!parse_param(str, object)) return parse_error();
      if (named_objects.find(object) == named_objects.end())
        return parse_error();
      auto& named_object = named_objects.at(object);
      for (auto& shape_id : named_object) {
        pbrt.shapes[shape_id].instances.push_back(
            flatten(ctx.stack.back().transform_start));
        pbrt.shapes[shape_id].instaends.push_back(
            flatten(ctx.stack.back().transform_end));
      }
    } else if (cmd == "ActiveTransform") {
      auto name = ""s;
      if (!parse_command(str, name)) return parse_error();
      if (name == "StartTime") {
        ctx.stack.back().active_transform_start = true;
        ctx.stack.back().active_transform_end   = false;
      } else if (name == "EndTime") {
        ctx.stack.back().active_transform_start = false;
        ctx.stack.back().active_transform_end   = true;
      } else if (name == "All") {
        ctx.stack.back().active_transform_start = true;
        ctx.stack.back().active_transform_end   = true;
      } else {
        throw std::out_of_range{"invalid command"};
      }
    } else if (cmd == "Transform") {
      auto xf = array<float, 16>{};
      if (!parse_param(str, xf)) return parse_error();
      set_transform(ctx.stack.back(), mat_to_frame(unflatten(xf)));
    } else if (cmd == "ConcatTransform") {
      auto xf = array<float, 16>{};
      if (!parse_param(str, xf)) return parse_error();
      concat_transform(ctx.stack.back(), mat_to_frame(unflatten(xf)));
    } else if (cmd == "Scale") {
      auto v = array<float, 3>{0, 0, 0};
      if (!parse_param(str, v)) return parse_error();
      concat_transform(ctx.stack.back(), scaling_frame(v));
    } else if (cmd == "Translate") {
      auto v = array<float, 3>{0, 0, 0};
      if (!parse_param(str, v)) return parse_error();
      concat_transform(ctx.stack.back(), translation_frame(v));
    } else if (cmd == "Rotate") {
      auto v = array<float, 4>{0, 0, 0, 0};
      if (!parse_param(str, v)) return parse_error();
      concat_transform(
          ctx.stack.back(), rotation_frame(array<float, 3>{v[1], v[2], v[3]},
                                v[0] * (float)M_PI / 180));
    } else if (cmd == "LookAt") {
      auto from = array<float, 3>{0, 0, 0}, to = array<float, 3>{0, 0, 0},
           up = array<float, 3>{0, 0, 0};
      if (!parse_param(str, from)) return parse_error();
      if (!parse_param(str, to)) return parse_error();
      if (!parse_param(str, up)) return parse_error();
      auto frame = lookat_frame(from, to, up, true);
      concat_transform(ctx.stack.back(), inverse(frame));
    } else if (cmd == "ReverseOrientation") {
      ctx.stack.back().reverse = !ctx.stack.back().reverse;
    } else if (cmd == "CoordinateSystem") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      ctx.coordsys[name].transform_start = ctx.stack.back().transform_start;
      ctx.coordsys[name].transform_end   = ctx.stack.back().transform_end;
    } else if (cmd == "CoordSysTransform") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      if (ctx.coordsys.find(name) != ctx.coordsys.end()) {
        ctx.stack.back().transform_start =
            ctx.coordsys.at(name).transform_start;
        ctx.stack.back().transform_end = ctx.coordsys.at(name).transform_end;
      }
    } else if (cmd == "Integrator") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Sampler") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "PixelFilter") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Film") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      auto film = pbrt_film{};
      if (!convert_film(film, command, filename)) return parse_error();
      ctx.film_resolution = film.resolution;
    } else if (cmd == "Accelerator") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Camera") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto& camera  = pbrt.cameras.emplace_back();
      if (!convert_camera(camera, command, ctx.film_resolution, filename))
        return parse_error();
    } else if (cmd == "Texture") {
      auto command  = pbrt_command{};
      auto comptype = ""s;
      auto str_     = string{str};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_param(str, comptype)) return parse_error();
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (!convert_texture(
              named_textures[command.name], command, named_textures, filename))
        return parse_error();
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        command     = pbrt_command{};
      command.name = "__unnamed__material__" + std::to_string(material_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (command.type.empty()) {
        ctx.stack.back().material = {};
      } else {
        ctx.stack.back().material = {};
        if (!convert_material(ctx.stack.back().material, command, texture_map,
                named_materials, named_textures, filename))
          return parse_error();
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (value.name == "type") command.type = value.value1s;
      if (!convert_material(named_materials[command.name], command, texture_map,
              named_materials, named_textures, filename))
        return parse_error();
    } else if (cmd == "NamedMaterial") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      if (named_materials.find(name) == named_materials.end())
        return parse_error();
      ctx.stack.back().material = named_materials.at(name);
    } else if (cmd == "Shape") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame  = ctx.stack.back().transform_start;
      command.frend  = ctx.stack.back().transform_end;
      auto& shape    = pbrt.shapes.emplace_back();
      auto  alphamap = ""s;
      if (!convert_shape(shape, command, alphamap, named_textures, ply_dirname,
              ply_meshes, filename))
        return parse_error();
      auto matkey = "?!!!?" + ctx.stack.back().material.name + "?!!!?" +
                    ctx.stack.back().arealight.name + "?!!!?" + alphamap;
      if (material_map.find(matkey) == material_map.end()) {
        auto& material    = pbrt.materials.emplace_back();
        material          = ctx.stack.back().material;
        material.name     = "material" + std::to_string(pbrt.materials.size());
        material.emission = ctx.stack.back().arealight.emission;
        // material.alpha_tex = alphamap;
        material_map[matkey] = (int)pbrt.materials.size() - 1;
      }
      shape.material = material_map.at(matkey);
      if (!ctx.cur_object.empty()) {
        named_objects[ctx.cur_object].push_back((int)pbrt.shapes.size() - 1);
        shape.instanced = true;
      }
    } else if (cmd == "AreaLightSource") {
      static auto arealight_id = 0;
      auto        command      = pbrt_command{};
      command.name = "__unnamed__arealight__" + std::to_string(arealight_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (!convert_arealight(ctx.stack.back().arealight, command, filename))
        return parse_error();
    } else if (cmd == "LightSource") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (command.type == "infinite") {
        auto& environment = pbrt.environments.emplace_back();
        if (!convert_environment(environment, command, texture_map, filename))
          return parse_error();
      } else {
        auto& light = pbrt.lights.emplace_back();
        if (!convert_light(light, command, filename)) return parse_error();
      }
    } else if (cmd == "MakeNamedMedium") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (command.name == "type") command.type = value.value1s;
      auto medium                 = pbrt_medium{};
      named_mediums[command.name] = medium;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      if (!parse_param(str, interior)) return parse_error();
      if (!parse_param(str, exterior)) return parse_error();
      ctx.stack.back().interior = named_mediums.at(interior);
      ctx.stack.back().exterior = named_mediums.at(exterior);
    } else if (cmd == "Include") {
      auto includename = ""s;
      if (!parse_param(str, includename)) return parse_error();
      if (!load_pbrt(path_join(path_dirname(filename), includename), pbrt,
              error, ctx, material_map, texture_map, named_materials,
              named_textures, named_mediums, named_objects, ply_dirname,
              ply_meshes))
        return dependent_error();
    } else {
      return parse_error();
    }
  }
  return true;
}

// load pbrt
bool load_pbrt(
    const string& filename, pbrt_model& pbrt, string& error, bool ply_meshes) {
  auto ctx             = pbrt_context{};
  auto material_map    = unordered_map<string, int>{};
  auto texture_map     = unordered_map<string, int>{};
  auto named_materials = unordered_map<string, pbrt_material>{{"", {}}};
  auto named_mediums   = unordered_map<string, pbrt_medium>{{"", {}}};
  auto named_textures  = unordered_map<string, pbrt_texture>{{"", {}}};
  auto named_objects   = unordered_map<string, vector<int>>{};
  if (!load_pbrt(filename, pbrt, error, ctx, material_map, texture_map,
          named_materials, named_textures, named_mediums, named_objects,
          path_dirname(filename), ply_meshes))
    return false;
  pbrt.textures.resize(texture_map.size());
  for (auto& [path, texture_id] : texture_map) {
    pbrt.textures[texture_id].filename = path;
  }
  return true;
}

static void format_value(string& str, const pbrt_value& value) {
  static auto type_labels = unordered_map<pbrt_type, string>{
      {pbrt_type::real, "float"},
      {pbrt_type::integer, "integer"},
      {pbrt_type::boolean, "bool"},
      {pbrt_type::string, "string"},
      {pbrt_type::point, "point"},
      {pbrt_type::normal, "normal"},
      {pbrt_type::vector, "vector"},
      {pbrt_type::texture, "texture"},
      {pbrt_type::color, "rgb"},
      {pbrt_type::point2, "point2"},
      {pbrt_type::vector2, "vector2"},
      {pbrt_type::spectrum, "spectrum"},
  };

  auto format_vector = [](string& str, auto& values) {
    str += "[ ";
    for (auto& value : values) {
      str += " ";
      format_value(str, value);
    }
    str += " ]";
  };

  format_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case pbrt_type::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_value(str, value.value1f);
      }
      break;
    case pbrt_type::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_value(str, value.value1i);
      }
      break;
    case pbrt_type::boolean:
      format_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case pbrt_type::string:
    case pbrt_type::texture: format_values(str, "\"{}\"", value.value1s); break;
    case pbrt_type::point:
    case pbrt_type::vector:
    case pbrt_type::normal:
    case pbrt_type::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_values(str, "[ {} ]", value.value3f);
      }
      break;
    case pbrt_type::spectrum: format_vector(str, value.vector1f); break;
    case pbrt_type::point2:
    case pbrt_type::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

static void format_value(string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

bool save_pbrt(const string& filename, const pbrt_model& pbrt, string& error,
    bool ply_meshes) {
  // buffer
  auto buffer = string{};

  // save comments
  format_values(buffer, "#\n");
  format_values(buffer, "# Written by Yocto/GL\n");
  format_values(buffer, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(buffer, "#\n\n");
  for (auto& comment : pbrt.comments) {
    format_values(buffer, "# {}\n", comment);
  }
  format_values(buffer, "\n");

  for (auto& camera : pbrt.cameras) {
    auto command = pbrt_command{};
    command.type = "image";
    command.values.push_back(
        make_pbrt_value("xresolution", camera.resolution[0]));
    command.values.push_back(
        make_pbrt_value("yresolution", camera.resolution[1]));
    command.values.push_back(make_pbrt_value("filename", "image.exr"s));
    format_values(buffer, "Film \"{}\" {}\n", command.type, command.values);
  }

  for (auto& camera : pbrt.cameras) {
    auto command  = pbrt_command{};
    command.type  = "perspective";
    auto cframe   = unflatten(camera.frame);
    command.frame = unflatten(camera.frame);
    command.values.push_back(make_pbrt_value(
        "fov", 2 * std::tan(0.036f / (2 * camera.lens)) * 180 / (float)M_PI));
    format_values(buffer, "LookAt {} {} {}\n", cframe[3],
        sub(cframe[3], cframe[2]), cframe[1]);
    format_values(buffer, "Camera \"{}\" {}\n", command.type, command.values);
  }

  format_values(buffer, "\nWorldBegin\n\n");

  for (auto& light : pbrt.lights) {
    auto command  = pbrt_command{};
    command.frame = unflatten(light.frame);
    if (light.distant) {
      command.type = "distance";
      command.values.push_back(make_pbrt_value("L", light.emission));
    } else {
      command.type = "point";
      command.values.push_back(make_pbrt_value("I", light.emission));
    }
    format_values(buffer, "AttributeBegin\n");
    format_values(buffer, "Transform {}\n", frame_to_mat(command.frame));
    format_values(
        buffer, "LightSource \"{}\" {}\n", command.type, command.values);
    format_values(buffer, "AttributeEnd\n");
  }

  for (auto& environment : pbrt.environments) {
    auto command  = pbrt_command{};
    command.frame = unflatten(environment.frame);
    command.type  = "infinite";
    command.values.push_back(make_pbrt_value("L", environment.emission));
    command.values.push_back(
        make_pbrt_value("mapname", environment.emission_tex));
    format_values(buffer, "AttributeBegin\n");
    format_values(buffer, "Transform {}\n", frame_to_mat(command.frame));
    format_values(
        buffer, "LightSource \"{}\" {}\n", command.type, command.values);
    format_values(buffer, "AttributeEnd\n");
  }

  auto reflectivity_to_eta1 = [](float reflectivity) {
    return (1 + std::sqrt(reflectivity)) / (1 - std::sqrt(reflectivity));
  };
  auto reflectivity_to_eta =
      [&](const array<float, 3>& reflectivity) -> array<float, 3> {
    return {reflectivity_to_eta1(reflectivity[0]),
        reflectivity_to_eta1(reflectivity[1]),
        reflectivity_to_eta1(reflectivity[2])};
  };

  for (auto& material : pbrt.materials) {
    auto command = pbrt_command{};
    switch (material.type) {
      case pbrt_mtype::matte: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
      } break;
      case pbrt_mtype::plastic: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
        command.values.push_back(
            make_pbrt_value("Ks", array<float, 3>{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", std::pow(material.roughness, 2.0f)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::metal: {
        command.type = "metal";
        command.values.push_back(
            make_pbrt_value("Kr", array<float, 3>{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", std::pow(material.roughness, 2.0f)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::thinglass: {
        command.type = "uber";
        command.values.push_back(
            make_pbrt_value("Ks", array<float, 3>{1, 1, 1}));
        command.values.push_back(make_pbrt_value("Kt", material.color));
        command.values.push_back(
            make_pbrt_value("roughness", std::pow(material.roughness, 2.0f)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::glass: {
        command.type = "glass";
        command.values.push_back(
            make_pbrt_value("Kr", array<float, 3>{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("Kt", array<float, 3>{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", std::pow(material.roughness, 2.0f)));
        command.values.push_back(make_pbrt_value("eta", material.ior));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::subsurface: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
      } break;
    }

    format_values(buffer,
        "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n", material.name,
        command.type, command.values);
  }

  auto object_id = 0;
  for (auto& shape : pbrt.shapes) {
    auto& material = pbrt.materials.at(shape.material);
    auto  command  = pbrt_command{};
    command.frame  = unflatten(shape.frame);
    if (ply_meshes) {
      command.type = "plymesh";
      command.values.push_back(make_pbrt_value("filename", shape.filename_));
    } else {
      command.type = "trianglemesh";
      command.values.push_back(make_pbrt_value("indices", shape.triangles));
      command.values.push_back(
          make_pbrt_value("P", shape.positions, pbrt_type::point));
      if (!shape.normals.empty())
        command.values.push_back(
            make_pbrt_value("N", shape.triangles, pbrt_type::normal));
      if (!shape.texcoords.empty())
        command.values.push_back(make_pbrt_value("uv", shape.texcoords));
    }
    if (ply_meshes) {
      auto ply = ply_model{};
      add_positions(ply, shape.positions);
      add_normals(ply, shape.normals);
      add_texcoords(ply, shape.texcoords);
      add_triangles(ply, shape.triangles);
      if (!save_ply(path_dirname(filename) + "/" + shape.filename_, ply, error))
        return false;
    }
    auto object = "object" + std::to_string(object_id++);
    if (!shape.instances.empty())
      format_values(buffer, "ObjectBegin \"{}\"\n", object);
    format_values(buffer, "AttributeBegin\n");
    format_values(
        buffer, "Transform {}\n", frame_to_mat(unflatten(shape.frame)));
    if (material.emission != array<float, 3>{0, 0, 0}) {
      auto acommand = pbrt_command{};
      acommand.type = "diffuse";
      acommand.values.push_back(make_pbrt_value("L", material.emission));
      format_values(buffer, "AreaLightSource \"{}\" {}\n", acommand.type,
          acommand.values);
    }
    format_values(buffer, "NamedMaterial \"{}\"\n", material.name);
    format_values(buffer, "Shape \"{}\" {}\n", command.type, command.values);
    format_values(buffer, "AttributeEnd\n");
    if (!shape.instances.empty()) format_values(buffer, "ObjectEnd\n");
    for (auto& iframe : shape.instances) {
      format_values(buffer, "AttributeBegin\n");
      format_values(buffer, "Transform {}\n", frame_to_mat(unflatten(iframe)));
      format_values(buffer, "ObjectInstance \"{}\"\n", object);
      format_values(buffer, "AttributeEnd\n");
    }
  }

  format_values(buffer, "\nWorldEnd\n\n");

  // save file
  if (!save_text(filename, buffer, error)) return false;

  // done
  return true;
}

}  // namespace yocto
