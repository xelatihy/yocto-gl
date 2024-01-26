//
// # Yocto/Animation: Animation operations
//
// Yocto/Animation provides basic utilities for working with animation.
//

//
// LICENSE:
//
// Copyright (c) 2024 -- 2024 Fabio Pellacini
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

#ifndef _YOCTO_ANIMATION_H_
#define _YOCTO_ANIMATION_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <utility>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternions to represent rotations
struct quat4f {
  vec3f v = {0, 0, 0};
  float w = 1;

  constexpr quat4f() : v{0, 0, 0}, w{1} {}
  constexpr quat4f(vec3f v_, float w_) : v{v_}, w{w_} {}
};

// Constants
constexpr auto identity_quat4f = quat4f{{0, 0, 0}, 1};

// Conversion to/from quaternions
inline quat4f             axisangle_to_quat(vec3f axis, float angle);
inline pair<vec3f, float> quat_to_axisangle(const quat4f& quaternion);

// Quaternion comparisons.
inline bool operator==(const quat4f& a, const quat4f& b);
inline bool operator!=(const quat4f& a, const quat4f& b);

// Quaternion operations
inline quat4f operator-(const quat4f& a);
inline quat4f operator+(const quat4f& a, const quat4f& b);
inline quat4f operator*(const quat4f& a, float b);
inline quat4f operator*(float a, const quat4f& b);
inline quat4f operator/(const quat4f& a, float b);
inline quat4f operator*(const quat4f& a, const quat4f& b);

// Quaternion operations
inline float  dot(const quat4f& a, const quat4f& b);
inline float  norm(const quat4f& a);
inline quat4f normalize(const quat4f& a);
inline float  length(const quat4f& a);
inline quat4f conjugate(const quat4f& a);
inline quat4f inverse(const quat4f& a);
inline float  uangle(const quat4f& a, const quat4f& b);
inline quat4f lerp(const quat4f& a, const quat4f& b, float t);
inline quat4f nlerp(const quat4f& a, const quat4f& b, float t);
inline quat4f slerp(const quat4f& a, const quat4f& b, float t);

}  // namespace yocto

// -----------------------------------------------------------------------------
// DUAL QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Dual quaternions
struct dualquat4f {
  quat4f q1 = {{0, 0, 0}, 1};
  quat4f q2 = {{0, 0, 0}, 0};

  constexpr dualquat4f() : q1{{0, 0, 0}, 1}, q2{{0, 0, 0}, 0} {}
  constexpr dualquat4f(const quat4f& q1_, const quat4f& q2_)
      : q1{q1_}, q2{q2_} {}
};

// Constants
constexpr auto identitydq4f = dualquat4f{{{0, 0, 0}, 1}, {{0, 0, 0}, 0}};

// Quaternion operations
inline dualquat4f operator+(const dualquat4f& a, const dualquat4f& b);
inline dualquat4f operator*(float a, const dualquat4f& b);
inline dualquat4f operator*(const dualquat4f& a, float b);
inline dualquat4f operator*(float a, const dualquat4f& b);
inline dualquat4f operator/(const dualquat4f& a, float b);
inline dualquat4f operator*(const dualquat4f& a, const dualquat4f& b);

// Quaternion operations
inline dualquat4f conjugate(const dualquat4f& a);
inline vec2f      norm(const dualquat4f& a);
inline dualquat4f normalize(const dualquat4f& a);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternion transformations
inline vec3f transform_point(const quat4f& a, vec3f b);
inline vec3f transform_vector(const quat4f& a, vec3f b);
inline vec3f transform_direction(const quat4f& a, vec3f b);

// Dual Quaternion transformations
inline vec3f transform_point(const dualquat4f& a, vec3f b);
inline vec3f transform_vector(const dualquat4f& a, vec3f b);
inline vec3f transform_direction(const dualquat4f& a, vec3f b);

// Quaternion transformations
inline frame3f rotation_frame(const quat4f& quat);

// Dual Quaternion transformations
inline dualquat4f translation_dualquat(vec3f translation);
inline dualquat4f rotation_dualquat(vec3f axis, float angle);
inline dualquat4f rotation_dualquat(const quat4f& rotation);
inline dualquat4f rigidtransform_dualquat(
    vec3f translation, const quat4f& rotation);
inline frame3f rigidtransform_frame(const dualquat4f& transform);

}  // namespace yocto

// -----------------------------------------------------------------------------
// INTERPOLATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience axis-angle rotation
using axis_angle3f = pair<vec3f, float>;

// Linear interpolation of translation, scaling and rotations.
// Rotations are interpolated using the exponential map.
inline frame3f interpolate_translation(vec3f a, vec3f b, float t);
inline frame3f interpolate_rotation(
    const axis_angle3f& a, const axis_angle3f& b, float t);
inline frame3f interpolate_rotation(const quat4f& a, const quat4f& b, float t);
inline frame3f interpolate_scaling(vec3f a, vec3f b, float t);

// Exponential map for matrices and quaternions
inline mat3f        exp_map(const axis_angle3f& rotation);
inline axis_angle3f log_map(const mat3f& rotation);
inline quat4f       exp_map(const quat4f& rotation);
inline quat4f       log_map(const quat4f& rotation);

}  // namespace yocto

// -----------------------------------------------------------------------------
// DEFORMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply linear skinning to vertex position and normals.
inline pair<vector<vec3f>, vector<vec3f>> skin_linear(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms);
// Update skinning
inline void skin_linear(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms);

// Apply dual quaternion skinning to vertex position and normals.
inline pair<vector<vec3f>, vector<vec3f>> skin_dualquat(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<dualquat4f>& xforms);
// Update skinning
inline void skin_dualquat(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<dualquat4f>& xforms);

// Apply glTF matrix skinning to vertex position and normals.
inline pair<vector<vec3f>, vector<vec3f>> skin_matrices(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms);
inline void skin_matrices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion to/from quaternions
inline quat4f axisangle_to_quat(vec3f axis, float angle) {
  return {sin(angle / 2) * normalize(axis), cos(angle / 2)};
}
inline pair<vec3f, float> quat_to_axisangle(const quat4f& quaternion) {
  return {
      normalize(quaternion.v), 2 * atan(length(quaternion.v), quaternion.w)};
}

// Quaternion comparisons.
inline bool operator==(const quat4f& a, const quat4f& b) {
  return a.v == b.v && a.w == b.w;
}
inline bool operator!=(const quat4f& a, const quat4f& b) {
  return a.v != b.v || a.w != b.w;
}

// Quaternion operatons
inline quat4f operator-(const quat4f& a) { return {-a.v, -a.w}; }
inline quat4f operator+(const quat4f& a, const quat4f& b) {
  return {a.v + b.v, a.w + b.w};
}
inline quat4f operator*(const quat4f& a, float b) { return {a.v * b, a.w * b}; }
inline quat4f operator*(float a, const quat4f& b) { return {a * b.v, a * b.w}; }
inline quat4f operator/(const quat4f& a, float b) { return {a.v / b, a.w / b}; }
inline quat4f operator*(const quat4f& a, const quat4f& b) {
  return {a.v * b.w + a.w * b.v + cross(a.v, b.v), a.w * b.w - dot(a.v, b.v)};
}

// Quaternion operations
inline float dot(const quat4f& a, const quat4f& b) {
  return dot(a.v, b.v) + a.w * b.w;
}
inline float  norm(const quat4f& a) { return sqrt(dot(a, a)); }
inline quat4f normalize(const quat4f& a) {
  auto n = norm(a);
  return (n != 0) ? a / n : a;
}
inline float  length(const quat4f& a) { return sqrt(dot(a, a)); }
inline quat4f conjugate(const quat4f& a) { return {-a.v, a.w}; }
inline quat4f inverse(const quat4f& a) { return conjugate(a) / dot(a, a); }
inline float  uangle(const quat4f& a, const quat4f& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
inline quat4f lerp(const quat4f& a, const quat4f& b, float t) {
  return a * (1 - t) + b * t;
}
inline quat4f nlerp(const quat4f& a, const quat4f& b, float t) {
  return normalize(lerp(a, b, t));
}
inline quat4f slerp(const quat4f& a, const quat4f& b, float t) {
  auto th = uangle(a, b);
  return th == 0
             ? a
             : a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
}

// Quaternion functions
inline quat4f exp(const quat4f& a) {
  return quat4f{normalize(a.v) * sin(length(a.v)), cos(length(a.v))} * exp(a.w);
}
inline quat4f log(const quat4f& a) {
  return quat4f{normalize(a.v) * acos(a.w / length(a)), log(length(a))};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DUAL QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternion operations
inline dualquat4f operator+(const dualquat4f& a, const dualquat4f& b) {
  return {a.q1 + b.q1, a.q2 + b.q2};
}
inline dualquat4f operator*(const dualquat4f& a, float b) {
  return {a.q1 * b, a.q2 * b};
}
inline dualquat4f operator*(float a, const dualquat4f& b) {
  return {a * b.q1, a * b.q2};
}
inline dualquat4f operator/(const dualquat4f& a, float b) {
  return {a.q1 / b, a.q2 / b};
}
inline dualquat4f operator*(const dualquat4f& a, const dualquat4f& b) {
  return {a.q1 * b.q1, a.q1 * b.q2 + a.q2 * b.q1};
}

// Quaternion operations
inline dualquat4f conjugate(const dualquat4f& a) {
  return {conjugate(a.q1), conjugate(a.q2)};
}
inline vec2f norm(const dualquat4f& a) {
  return {norm(a.q1), dot(a.q1, a.q2) / length(a.q1)};
}
inline dualquat4f normalize(const dualquat4f& a) {
  auto l = norm(a.q1);
  return {a.q1 / l, a.q2 / l};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternion transformations
inline vec3f transform_point(const dualquat4f& a, vec3f b) {
  return (a * dualquat4f{{zero3f, 1}, {b, 0}} *
          dualquat4f{conjugate(a.q1), -conjugate(a.q2)})
      .q2.v;
}
inline vec3f transform_vector(const dualquat4f& a, vec3f b) {
  return (a * dualquat4f{{zero3f, 0}, {b, 0}} *
          dualquat4f{conjugate(a.q1), -conjugate(a.q2)})
      .q2.v;
}
inline vec3f transform_direction(const dualquat4f& a, vec3f b) {
  return normalize(transform_vector(a, b));
}

// Quaternion transformations
inline dualquat4f translation_dualquat(vec3f translation) {
  return {{zero3f, 1}, {translation / 2, 0}};
}
inline dualquat4f rotation_dualquat(vec3f axis, float angle) {
  return {{normalize(axis) * sin(angle / 2), cos(angle / 2)}, {zero3f, 0}};
}
inline dualquat4f rotation_dualquat(const quat4f& rotation) {
  return {rotation, {zero3f, 0}};
}
inline dualquat4f rigidtransform_dualquat(
    vec3f translation, const quat4f& rotation) {
  return translation_dualquat(translation) * rotation_dualquat(rotation);
}
inline frame3f rigidtransform_frame(const dualquat4f& transform) {
  return translation_frame(transform.q2.v) * rotation_frame(transform.q1);
}

inline frame3f rotation_frame(const quat4f& quat) {
  auto v = vec4f{quat.v.x, quat.v.y, quat.v.z, quat.w};  // convenience
  return {{v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
              (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
      {(v.x * v.y - v.z * v.w) * 2,
          v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
          (v.y * v.z + v.x * v.w) * 2},
      {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
          v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z},
      {0, 0, 0}};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INTERPOLATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Cross product matrix
inline mat3f cross_prod_mat(vec3f a) {
  return {{0, a.z, -a.y}, {-a.z, 0, a.x}, {a.y, -a.x, 0}};
}

// Exponential map for matrices
inline mat3f exp_map(const axis_angle3f& rotation) {
  auto [axis, angle] = rotation;
  auto cross_mat     = cross_prod_mat(axis);
  return identity3x3f + cross_mat * sin(angle) +
         (cross_mat * cross_mat) * (1 - cos(angle));
}
inline axis_angle3f log_map(const mat3f& rotation) {
  auto angle = acos(clamp((trace(rotation) - 1) / 2, (float)-1, (float)1));
  if (angle == 0) return {{0, 0, 1}, 0};
  return {vec3f{rotation.z.y - rotation.y.z, rotation.x.z - rotation.z.x,
              rotation.y.x - rotation.x.y} /
              (2 * sin(angle)),
      angle};
}
inline quat4f exp_map(const quat4f& q) {
  return quat4f{normalize(q.v) * sin(length(q.v)), cos(length(q.v))} * exp(q.w);
}
inline quat4f log_map(const quat4f& q) {
  return quat4f{normalize(q.v) * acos(q.w / length(q)), log(length(q))};
}

// Linear interpolation of translation, scaling and rotations.
// Rotations are interpolated using the exponential map.
inline frame3f interpolate_translation(vec3f a, vec3f b, float t) {
  return translation_frame(lerp(a, b, t));
}
inline frame3f interpolate_scaling(vec3f a, vec3f b, float t) {
  return scaling_frame(lerp(a, b, t));
}
inline frame3f interpolate_rotation(
    const axis_angle3f& a, const axis_angle3f& b, float t) {
  auto ma = rotation(rotation_frame(a)), mb = rotation(rotation_frame(b));
  auto [axis, angle] = log_map(mb * transpose(ma));
  // TODO: fix the sign of the axis
  auto rotation = exp_map(axis_angle3f{-axis, angle * t}) * ma;
  return {rotation, zero3f};
}
inline frame3f interpolate_rotation(const quat4f& a, const quat4f& b, float t) {
  return rotation_frame(exp_map(log_map(b * conjugate(a)) * t) * a);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DEFORMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply skinning
inline pair<vector<vec3f>, vector<vec3f>> skin_linear(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i : range(positions.size())) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
  }
  for (auto i : range(normals.size())) {
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning
inline void skin_linear(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i : range(positions.size())) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
  }
  for (auto i : range(normals.size())) {
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
}

// Apply skinning
inline pair<vector<vec3f>, vector<vec3f>> skin_dualquat(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<dualquat4f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i : range(positions.size())) {
    auto dualquat        = normalize(xforms[joints[i].x] * weights[i].x +
                                     xforms[joints[i].y] * weights[i].y +
                                     xforms[joints[i].z] * weights[i].z +
                                     xforms[joints[i].w] * weights[i].w);
    skinned_positions[i] = transform_point(dualquat, positions[i]);
    if (!normals.empty())
      skinned_normals[i] = transform_direction(dualquat, normals[i]);
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning
inline void skin_dualquat(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<dualquat4f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i : range(positions.size())) {
    auto dualquat        = normalize(xforms[joints[i].x] * weights[i].x +
                                     xforms[joints[i].y] * weights[i].y +
                                     xforms[joints[i].z] * weights[i].z +
                                     xforms[joints[i].w] * weights[i].w);
    skinned_positions[i] = transform_point(dualquat, positions[i]);
    if (!normals.empty())
      skinned_normals[i] = transform_direction(dualquat, normals[i]);
  }
}

// Apply skinning as specified in Khronos glTF
inline pair<vector<vec3f>, vector<vec3f>> skin_matrices(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i : range(positions.size())) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning as specified in Khronos glTF
inline void skin_matrices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i : range(positions.size())) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Find the first keyframe value that is greater than the argument.
inline int keyframe_index(const vector<float>& times, const float& time) {
  for (auto i : range(times.size()))
    if (times[i] > time) return (int)i;
  return (int)times.size();
}

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T keyframe_step(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f keyframe_slerp(
    const vector<float>& times, const vector<vec4f>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline T keyframe_linear(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evaluates a keyframed value using Bezier interpolation.
template <typename T>
inline T keyframe_bezier(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return interpolate_bezier(
      vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace yocto

#endif
