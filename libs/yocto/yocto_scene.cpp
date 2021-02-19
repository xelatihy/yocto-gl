//
// Implementation for Yocto/Scene Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "ext/stb_image_resize.h"
#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_noise.h"
#include "yocto_parallel.h"
#include "yocto_sceneio.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// image creation
image_data make_image(int width, int height, bool linear, bool as_byte) {
  if (!as_byte) {
    return image_data{width, height, linear,
        vector<vec4f>(width * height, vec4f{0, 0, 0, 0}), {}};
  } else {
    return image_data{width, height, linear, {},
        vector<vec4b>(width * height, vec4b{0, 0, 0, 0})};
  }
}
image_data make_image(int width, int height, bool linear, const vec4f* data) {
  return image_data{
      width, height, linear, vector<vec4f>(data, data + width * height), {}};
}
image_data make_image(int width, int height, bool linear, const vec4b* data) {
  return image_data{
      width, height, linear, {}, vector<vec4b>(data, data + width * height)};
}

// equality
bool operator==(const image_data& a, const image_data& b) {
  return a.width == b.width && a.height == b.height && a.linear == b.linear &&
         a.pixelsf == b.pixelsf && a.pixelsb == b.pixelsb;
}
bool operator!=(const image_data& a, const image_data& b) {
  return a.width != b.width || a.height != b.height || a.linear != b.linear ||
         a.pixelsf != b.pixelsf || a.pixelsb != b.pixelsb;
}

// swap
void swap(image_data& a, image_data& b) {
  std::swap(a.width, b.width);
  std::swap(a.height, b.height);
  std::swap(a.linear, b.linear);
  std::swap(a.pixelsf, b.pixelsf);
  std::swap(a.pixelsb, b.pixelsb);
}

// pixel access
vec4f get_pixel(const image_data& image, int i, int j) {
  if (!image.pixelsf.empty()) {
    return image.pixelsf[j * image.width + i];
  } else {
    return byte_to_float(image.pixelsb[j * image.width + i]);
  }
}
void set_pixel(image_data& image, int i, int j, const vec4f& pixel) {
  if (!image.pixelsf.empty()) {
    image.pixelsf[j * image.width + i] = pixel;
  } else {
    image.pixelsb[j * image.width + i] = float_to_byte(pixel);
  }
}

// conversions
image_data convert_image(const image_data& image, bool linear, bool as_byte) {
  if (image.linear == linear && !image.pixelsb.empty() == as_byte) return image;
  auto result = make_image(image.width, image.height, linear, as_byte);
  convert_image(result, image);
  return result;
}
void convert_image(image_data& result, const image_data& image) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image have to be the same size"};
  if (image.linear == result.linear &&
      image.pixelsf.empty() == result.pixelsf.empty()) {
    result.pixelsf = image.pixelsf;
    result.pixelsb = image.pixelsb;
  } else if (image.linear == result.linear) {
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        set_pixel(result, i, j, get_pixel(image, i, j));
      }
    }
  } else {
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        auto color     = get_pixel(image, i, j);
        auto converted = image.linear ? rgb_to_srgb(color) : srgb_to_rgb(color);
        set_pixel(result, i, j, converted);
      }
    }
  }
}

// Evaluates an image at a point `uv`.
vec4f eval_image(const image_data& image, const vec2f& uv, bool as_linear,
    bool no_interpolation, bool clamp_to_edge) {
  if (image.width == 0 || image.height == 0) return {0, 0, 0, 0};

  // get image width/height
  auto size = vec2i{image.width, image.height};

  // get coordinates normalized for tiling
  auto s = 0.0f, t = 0.0f;
  if (clamp_to_edge) {
    s = clamp(uv.x, 0.0f, 1.0f) * size.x;
    t = clamp(uv.y, 0.0f, 1.0f) * size.y;
  } else {
    s = fmod(uv.x, 1.0f) * size.x;
    if (s < 0) s += size.x;
    t = fmod(uv.y, 1.0f) * size.y;
    if (t < 0) t += size.y;
  }

  // get image coordinates and residuals
  auto i = clamp((int)s, 0, size.x - 1), j = clamp((int)t, 0, size.y - 1);
  auto ii = (i + 1) % size.x, jj = (j + 1) % size.y;
  auto u = s - i, v = t - j;

  if (no_interpolation) {
    if (as_linear && !image.linear) {
      return srgb_to_rgb(get_pixel(image, i, j));
    } else {
      return get_pixel(image, i, j);
    }
  } else {
    // handle interpolation
    if (as_linear && !image.linear) {
      return srgb_to_rgb(get_pixel(image, i, j)) * (1 - u) * (1 - v) +
             srgb_to_rgb(get_pixel(image, i, jj)) * (1 - u) * v +
             srgb_to_rgb(get_pixel(image, ii, j)) * u * (1 - v) +
             srgb_to_rgb(get_pixel(image, ii, jj)) * u * v;
    } else {
      return get_pixel(image, i, j) * (1 - u) * (1 - v) +
             get_pixel(image, i, jj) * (1 - u) * v +
             get_pixel(image, ii, j) * u * (1 - v) +
             get_pixel(image, ii, jj) * u * v;
    }
  }
}

// Apply tone mapping returning a float or byte image.
image_data tonemap_image(
    const image_data& image, float exposure, bool filmic, bool as_byte) {
  if (!image.linear) return image;
  auto result = make_image(image.width, image.height, false, as_byte);
  for (auto idx = 0; idx < image.width * image.height; idx++) {
    result.pixelsb[idx] = float_to_byte(
        tonemap(image.pixelsf[idx], exposure, filmic, true));
  }
  return result;
}

// Apply tone mapping. If the input image is an ldr, does nothing.
void tonemap_image(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (!image.linear) throw std::invalid_argument{"hdr expected"};
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto hdr = get_pixel(image, i, j);
      auto ldr = tonemap(hdr, exposure, filmic);
      set_pixel(result, i, j, ldr);
    }
  }
}
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (!image.linear) throw std::invalid_argument{"hdr expected"};
  parallel_for(image.width, image.height,
      [&result, &image, exposure, filmic](int i, int j) {
        auto hdr = get_pixel(image, i, j);
        auto ldr = tonemap(hdr, exposure, filmic);
        set_pixel(result, i, j, ldr);
      });
}

// Resize an image.
image_data resize_image(
    const image_data& image, int res_width, int res_height) {
  if (res_width == 0 && res_height == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (res_height == 0) {
    res_height = (int)round(
        res_width * (double)image.height / (double)image.width);
  } else if (res_width == 0) {
    res_width = (int)round(
        res_height * (double)image.width / (double)image.height);
  }
  if (!image.pixelsf.empty()) {
    auto result = make_image(res_width, res_height, image.linear, true);
    stbir_resize_uint8_generic((byte*)image.pixelsb.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4b) * image.width),
        (byte*)result.pixelsb.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4b) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  } else {
    auto result = make_image(res_width, res_height, image.linear, false);
    stbir_resize_float_generic((float*)image.pixelsf.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4f) * image.width),
        (float*)result.pixelsf.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4f) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  }
}

// Compute the difference between two images.
image_data image_difference(
    const image_data& image1, const image_data& image2, bool display) {
  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    throw std::invalid_argument{"image sizes are different"};
  }

  // check types
  if (!image1.pixelsf.empty() != !image2.pixelsf.empty() ||
      !image1.pixelsf.empty() != !image2.pixelsf.empty()) {
    throw std::invalid_argument{"image types are different"};
  }

  // check types
  if (image1.linear != image2.linear || !image1.linear != !image2.linear) {
    throw std::invalid_argument{"image types are different"};
  }

  // compute diff
  auto difference = make_image(
      image1.width, image1.height, image1.linear, false);
  for (auto j = 0; j < image1.height; j++) {
    for (auto i = 0; i < image1.width; i++) {
      auto diff = abs(get_pixel(image1, i, j) - get_pixel(image2, i, j));
      if (display) {
        auto d = max(diff);
        set_pixel(difference, i, j, {d, d, d, 1});
      } else {
        set_pixel(difference, i, j, diff);
      }
    }
  }
  return difference;
}

void set_region(image_data& image, const image_data& region, int x, int y) {
  for (auto j = 0; j < region.height; j++) {
    for (auto i = 0; i < region.width; i++) {
      set_pixel(image, i + x, j + y, get_pixel(region, i, j));
    }
  }
}

void get_region(image_data& region, const image_data& image, int x, int y,
    int width, int height) {
  if (region.width != width || region.height != height) {
    region = make_image(width, height, image.linear, !image.pixelsf.empty());
  }
  for (auto j = 0; j < height; j++) {
    for (auto i = 0; i < width; i++) {
      set_pixel(region, i, j, get_pixel(region, i + x, j + y));
    }
  }
}

// Apply color grading from a linear or srgb color to an srgb color.
vec4f colorgradeb(
    const vec4f& color, bool linear, const colorgrade_params& params) {
  auto rgb   = xyz(color);
  auto alpha = color.w;
  if (linear) {
    if (params.exposure != 0) rgb *= exp2(params.exposure);
    if (params.tint != vec3f{1, 1, 1}) rgb *= params.tint;
    if (params.lincontrast != 0.5f)
      rgb = lincontrast(rgb, params.lincontrast, 0.18f);
    if (params.logcontrast != 0.5f)
      rgb = logcontrast(rgb, params.logcontrast, 0.18f);
    if (params.linsaturation != 0.5f) rgb = saturate(rgb, params.linsaturation);
    if (params.filmic) rgb = tonemap_filmic(rgb);
    if (params.srgb) rgb = rgb_to_srgb(rgb);
  }
  if (params.contrast != 0.5f) rgb = contrast(rgb, params.contrast);
  if (params.saturation != 0.5f) rgb = saturate(rgb, params.saturation);
  if (params.shadows != 0.5f || params.midtones != 0.5f ||
      params.highlights != 0.5f || params.shadows_color != vec3f{1, 1, 1} ||
      params.midtones_color != vec3f{1, 1, 1} ||
      params.highlights_color != vec3f{1, 1, 1}) {
    auto lift  = params.shadows_color;
    auto gamma = params.midtones_color;
    auto gain  = params.highlights_color;
    lift       = lift - mean(lift) + params.shadows - (float)0.5;
    gain       = gain - mean(gain) + params.highlights + (float)0.5;
    auto grey  = gamma - mean(gamma) + params.midtones;
    gamma      = log(((float)0.5 - lift) / (gain - lift)) / log(grey);
    // apply_image
    auto lerp_value = clamp(pow(rgb, 1 / gamma), 0, 1);
    rgb             = gain * lerp_value + lift * (1 - lerp_value);
  }
  return vec4f{rgb.x, rgb.y, rgb.z, alpha};
}

// Color grade an hsr or ldr image to an ldr image.
image_data colorgrade_image(
    const image_data& image, const colorgrade_params& params, bool as_byte) {
  auto result = make_image(image.width, image.height, false, as_byte);
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto color  = get_pixel(image, i, j);
      auto graded = colorgrade(color, image.linear, params);
      set_pixel(result, i, j, graded);
    }
  }
  return result;
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image(image_data& result, const image_data& image,
    const colorgrade_params& params) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"non linear expected"};
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto color  = get_pixel(image, i, j);
      auto graded = colorgrade(color, image.linear, params);
      set_pixel(result, i, j, graded);
    }
  }
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(image_data& result, const image_data& image,
    const colorgrade_params& params) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"non linear expected"};
  parallel_for(
      image.width, image.height, [&result, &image, &params](int i, int j) {
        auto color  = get_pixel(image, i, j);
        auto graded = colorgrade(color, image.linear, params);
        set_pixel(result, i, j, graded);
      });
}

// determine white balance colors
vec4f compute_white_balance(const image_data& image) {
  auto rgb = vec3f{0, 0, 0};
  for (auto j = 0; image.height; j++) {
    for (auto i = 0; image.width; i++) {
      rgb += xyz(get_pixel(image, i, j));
    }
  }
  if (rgb == vec3f{0, 0, 0}) return {0, 0, 0, 1};
  rgb /= max(rgb);
  return {rgb.x, rgb.y, rgb.z, 1};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CAMERA PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera for yimg::image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera(
    const scene_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  if (!camera.orthographic) {
    auto q = vec3f{film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f),
        camera.lens};
    // ray direction through the lens center
    auto dc = -normalize(q);
    // point on the lens
    auto e = vec3f{
        lens_uv.x * camera.aperture / 2, lens_uv.y * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = dc * camera.focus / abs(dc.z);
    // correct ray direction to account for camera focusing
    auto d = normalize(p - e);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  } else {
    auto scale = 1 / camera.lens;
    auto q     = vec3f{film.x * (0.5f - image_uv.x) * scale,
        film.y * (image_uv.y - 0.5f) * scale, camera.lens};
    // point on the lens
    auto e = vec3f{-q.x, -q.y, 0} + vec3f{lens_uv.x * camera.aperture / 2,
                                        lens_uv.y * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = vec3f{-q.x, -q.y, -camera.focus};
    // correct ray direction to account for camera focusing
    auto d = normalize(p - e);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate a texture
vec4f eval_texture(const scene_texture& texture, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge) {
  return eval_image(texture, uv, as_linear, no_interpolation, clamp_to_edge);
}

// Helpers
vec4f eval_texture(const scene_scene& scene, texture_handle texture,
    const vec2f& uv, bool ldr_as_linear, bool no_interpolation,
    bool clamp_to_edge) {
  if (texture == invalid_handle) return {1, 1, 1, 1};
  return eval_texture(
      scene.textures[texture], uv, ldr_as_linear, no_interpolation);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATERIAL PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// constant values
static const auto min_roughness = 0.03f * 0.03f;

// Evaluate material
material_point eval_material(const scene_scene& scene,
    const scene_material& material, const vec2f& texcoord,
    const vec4f& color_shp) {
  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::glass ||
      material.type == material_type::volume ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::metallic ||
      point.type == material_type::plastic) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// check if a material is a delta or volumetric
bool is_delta(const scene_material& material) {
  return (material.type == material_type::metal && material.roughness == 0) ||
         (material.type == material_type::glass && material.roughness == 0) ||
         (material.type == material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == material_type::volume);
}
bool is_volumetric(const scene_material& material) {
  return material.type == material_type::glass ||
         material.type == material_type::volume ||
         material.type == material_type::subsurface;
}

// check if a brdf is a delta
bool is_delta(const material_point& material) {
  return (material.type == material_type::metal && material.roughness == 0) ||
         (material.type == material_type::glass && material.roughness == 0) ||
         (material.type == material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == material_type::volume);
}
bool has_volume(const material_point& material) {
  return material.type == material_type::glass ||
         material.type == material_type::volume ||
         material.type == material_type::subsurface;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const shape_data& shape, int element, const vec2f& uv) {
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.positions[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.positions[line.x], shape.positions[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return normalize(shape.normals[point]);
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return normalize(
        interpolate_line(shape.normals[line.x], shape.normals[line.y], uv.x));
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return normalize(interpolate_triangle(shape.normals[triangle.x],
        shape.normals[triangle.y], shape.normals[triangle.z], uv));
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec3f eval_tangent(const shape_data& shape, int element, const vec2f& uv) {
  return eval_normal(shape, element, uv);
}

vec2f eval_texcoord(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.texcoords[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.texcoords[line.x], shape.texcoords[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.texcoords[triangle.x],
        shape.texcoords[triangle.y], shape.texcoords[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

vec4f eval_color(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.colors[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.colors[line.x], shape.colors[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.colors[triangle.x],
        shape.colors[triangle.y], shape.colors[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.colors[quad.x], shape.colors[quad.y],
        shape.colors[quad.z], shape.colors[quad.w], uv);
  } else {
    return {0, 0};
  }
}

float eval_radius(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.radius.empty()) return 0;
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.radius[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.radius[line.x], shape.radius[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.radius[triangle.x],
        shape.radius[triangle.y], shape.radius[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.radius[quad.x], shape.radius[quad.y],
        shape.radius[quad.z], shape.radius[quad.w], uv);
  } else {
    return 0;
  }
}

// Evaluate element normals
vec3f eval_element_normal(const shape_data& shape, int element) {
  if (!shape.points.empty()) {
    return {0, 0, 1};
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return line_tangent(shape.positions[line.x], shape.positions[line.y]);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return triangle_normal(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z]);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const shape_data& shape) {
  if (!shape.points.empty()) {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    return lines_tangents(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return triangles_normals(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return quads_normals(shape.quads, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const shape_data& shape) {
  if (!shape.points.empty()) {
    normals.assign(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    lines_tangents(normals, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    triangles_normals(normals, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    quads_normals(normals, shape.quads, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Shape sampling
vector<float> sample_shape_cdf(const shape_data& shape) {
  if (!shape.points.empty()) {
    return sample_points_cdf((int)shape.points.size());
  } else if (!shape.lines.empty()) {
    return sample_lines_cdf(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return sample_triangles_cdf(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return sample_quads_cdf(shape.quads, shape.positions);
  } else {
    return sample_points_cdf((int)shape.positions.size());
  }
}

void sample_shape_cdf(vector<float>& cdf, const shape_data& shape) {
  if (!shape.points.empty()) {
    sample_points_cdf(cdf, (int)shape.points.size());
  } else if (!shape.lines.empty()) {
    sample_lines_cdf(cdf, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    sample_triangles_cdf(cdf, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    sample_quads_cdf(cdf, shape.quads, shape.positions);
  } else {
    sample_points_cdf(cdf, (int)shape.positions.size());
  }
}

shape_point sample_shape(const shape_data& shape, const vector<float>& cdf,
    float rn, const vec2f& ruv) {
  if (!shape.points.empty()) {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  } else if (!shape.lines.empty()) {
    auto [element, u] = sample_lines(cdf, rn, ruv.x);
    return {element, {u, 0}};
  } else if (!shape.triangles.empty()) {
    auto [element, uv] = sample_triangles(cdf, rn, ruv);
    return {element, uv};
  } else if (!shape.quads.empty()) {
    auto [element, uv] = sample_quads(cdf, rn, ruv);
    return {element, uv};
  } else {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  }
}

vector<shape_point> sample_shape(
    const shape_data& shape, int num_samples, uint64_t seed) {
  auto cdf    = sample_shape_cdf(shape);
  auto points = vector<shape_point>(num_samples);
  auto rng    = make_rng(seed);
  for (auto& point : points) {
    point = sample_shape(shape, cdf, rand1f(rng), rand2f(rng));
  }
  return points;
}

// Conversions
shape_data quads_to_triangles(const shape_data& shape) {
  auto result = shape;
  quads_to_triangles(result, result);
  return result;
}
void quads_to_triangles(shape_data& result, const shape_data& shape) {
  result.triangles = quads_to_triangles(shape.quads);
  result.quads     = {};
}

// Subdivision
shape_data subdivide_shape(
    const shape_data& shape, int subdivisions, bool catmullclark) {
  // This should probably be reimplemented in a faster fashion,
  // but how it is not obvious
  auto subdivided = shape_data{};
  if (!shape.points.empty()) {
    // nothing to do
  } else if (!shape.lines.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_lines(
        shape.lines, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_lines(
        shape.lines, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_lines(
        shape.lines, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_lines(
        shape.lines, shape.radius, subdivisions);
    std::tie(subdivided.lines, subdivided.positions) = subdivide_lines(
        shape.lines, shape.positions, subdivisions);
  } else if (!shape.triangles.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_triangles(
        shape.triangles, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_triangles(
        shape.triangles, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_triangles(
        shape.triangles, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_triangles(
        shape.triangles, shape.radius, subdivisions);
    std::tie(subdivided.triangles, subdivided.positions) = subdivide_triangles(
        shape.triangles, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && !catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_quads(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_quads(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_quads(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_quads(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_quads(
        shape.quads, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_catmullclark(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_catmullclark(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_catmullclark(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_catmullclark(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_catmullclark(
        shape.quads, shape.positions, subdivisions);
  } else {
    // empty shape
  }
  return subdivided;
}

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, const vec2f& uv) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadsnorm[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec2f eval_texcoord(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadstexcoord[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

// Evaluate element normals
vec3f eval_element_normal(const fvshape_data& shape, int element) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const fvshape_data& shape) {
  if (!shape.quadspos.empty()) {
    return quads_normals(shape.quadspos, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const fvshape_data& shape) {
  if (!shape.quadspos.empty()) {
    quads_normals(normals, shape.quadspos, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Conversions
shape_data fvshape_to_shape(const fvshape_data& fvshape, bool as_triangles) {
  auto shape = shape_data{};
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, fvshape.quadspos, fvshape.quadsnorm,
      fvshape.quadstexcoord, fvshape.positions, fvshape.normals,
      fvshape.texcoords);
  return shape;
}
fvshape_data shape_to_fvshape(const shape_data& shape) {
  if (!shape.points.empty() || !shape.lines.empty())
    throw std::invalid_argument{"cannor convert shape"};
  auto fvshape          = fvshape_data{};
  fvshape.positions     = shape.positions;
  fvshape.normals       = shape.normals;
  fvshape.texcoords     = shape.texcoords;
  fvshape.quadspos      = !shape.quads.empty() ? shape.quads
                                               : triangles_to_quads(shape.triangles);
  fvshape.quadsnorm     = !shape.normals.empty() ? fvshape.quadspos
                                                 : vector<vec4i>{};
  fvshape.quadstexcoord = !shape.texcoords.empty() ? fvshape.quadspos
                                                   : vector<vec4i>{};
  return fvshape;
}

// Subdivision
fvshape_data subdivide_fvshape(
    const fvshape_data& shape, int subdivisions, bool catmullclark) {
  auto subdivided = fvshape_data{};
  if (!catmullclark) {
    std::tie(subdivided.quadspos, subdivided.positions) = subdivide_quads(
        shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_quads(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) = subdivide_quads(
        shape.quadstexcoord, shape.texcoords, subdivisions);
  } else {
    std::tie(subdivided.quadspos, subdivided.positions) =
        subdivide_catmullclark(shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_catmullclark(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) =
        subdivide_catmullclark(
            shape.quadstexcoord, shape.texcoords, subdivisions, true);
  }
  return subdivided;
}

vector<string> shape_stats(const shape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("points:       " + format(shape.points.size()));
  stats.push_back("lines:        " + format(shape.lines.size()));
  stats.push_back("triangles:    " + format(shape.triangles.size()));
  stats.push_back("quads:        " + format(shape.quads.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("colors:       " + format(shape.colors.size()));
  stats.push_back("radius:       " + format(shape.radius.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

vector<string> fvshape_stats(const fvshape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("fvquads:      " + format(shape.quadspos.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INSTANCE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Eval position
vec3f eval_position(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_point(
        instance.frame, interpolate_triangle(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z], uv));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_point(instance.frame,
        interpolate_quad(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w], uv));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_point(instance.frame,
        interpolate_line(shape.positions[l.x], shape.positions[l.y], uv.x));
  } else if (!shape.points.empty()) {
    return transform_point(
        instance.frame, shape.positions[shape.points[element]]);
  } else {
    return zero3f;
  }
}

// Shape element normal.
vec3f eval_element_normal(
    const scene_scene& scene, const scene_instance& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, triangle_normal(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z]));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_normal(
        instance.frame, quad_normal(shape.positions[q.x], shape.positions[q.y],
                            shape.positions[q.z], shape.positions[q.w]));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_normal(instance.frame,
        line_tangent(shape.positions[l.x], shape.positions[l.y]));
  } else if (!shape.points.empty()) {
    return {0, 0, 1};
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
vec3f eval_normal(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.normals.empty())
    return eval_element_normal(scene, instance, element);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, normalize(interpolate_triangle(shape.normals[t.x],
                            shape.normals[t.y], shape.normals[t.z], uv)));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_normal(instance.frame,
        normalize(interpolate_quad(shape.normals[q.x], shape.normals[q.y],
            shape.normals[q.z], shape.normals[q.w], uv)));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_normal(instance.frame,
        normalize(
            interpolate_line(shape.normals[l.x], shape.normals[l.y], uv.x)));
  } else if (!shape.points.empty()) {
    return transform_normal(
        instance.frame, normalize(shape.normals[shape.points[element]]));
  } else {
    return zero3f;
  }
}

// Eval texcoord
vec2f eval_texcoord(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.texcoords.empty()) return uv;
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.texcoords[t.x], shape.texcoords[t.y], shape.texcoords[t.z], uv);
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return interpolate_quad(shape.texcoords[q.x], shape.texcoords[q.y],
        shape.texcoords[q.z], shape.texcoords[q.w], uv);
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return interpolate_line(shape.texcoords[l.x], shape.texcoords[l.y], uv.x);
  } else if (!shape.points.empty()) {
    return shape.texcoords[shape.points[element]];
  } else {
    return zero2f;
  }
}

#if 0
// Shape element normal.
static pair<vec3f, vec3f> eval_tangents(
    const trace_shape& shape, int element, const vec2f& uv) {
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    if (shape.texcoords.empty()) {
      return triangle_tangents_fromuv(shape.positions[t.x],
          shape.positions[t.y], shape.positions[t.z], {0, 0}, {1, 0}, {0, 1});
    } else {
      return triangle_tangents_fromuv(shape.positions[t.x],
          shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
          shape.texcoords[t.y], shape.texcoords[t.z]);
    }
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    if (shape.texcoords.empty()) {
      return quad_tangents_fromuv(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w], {0, 0}, {1, 0}, {0, 1},
          {1, 1}, uv);
    } else {
      return quad_tangents_fromuv(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w], shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], shape.texcoords[q.w],
          uv);
    }
  } else {
    return {zero3f, zero3f};
  }
}
#endif

// Shape element normal.
pair<vec3f, vec3f> eval_element_tangents(
    const scene_scene& scene, const scene_instance& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty() && !shape.texcoords.empty()) {
    auto t        = shape.triangles[element];
    auto [tu, tv] = triangle_tangents_fromuv(shape.positions[t.x],
        shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
        shape.texcoords[t.y], shape.texcoords[t.z]);
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else if (!shape.quads.empty() && !shape.texcoords.empty()) {
    auto q        = shape.quads[element];
    auto [tu, tv] = quad_tangents_fromuv(shape.positions[q.x],
        shape.positions[q.y], shape.positions[q.z], shape.positions[q.w],
        shape.texcoords[q.x], shape.texcoords[q.y], shape.texcoords[q.z],
        shape.texcoords[q.w], {0, 0});
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else {
    return {};
  }
}

vec3f eval_normalmap(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalid_handle &&
      (!shape.triangles.empty() || !shape.quads.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, false));
    auto [tu, tv]    = eval_element_tangents(scene, instance, element);
    auto frame       = frame3f{tu, tv, normal, zero3f};
    frame.x          = orthonormalize(frame.x, frame.z);
    frame.y          = normalize(cross(frame.z, frame.x));
    auto flip_v      = dot(frame.y, tv) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading normal
vec3f eval_shading_normal(const scene_scene& scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalid_handle) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::glass) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else if (!shape.lines.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    return orthonormalize(outgoing, normal);
  } else if (!shape.points.empty()) {
    return outgoing;
  } else {
    return zero3f;
  }
}

// Eval color
vec4f eval_color(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.colors[t.x], shape.colors[t.y], shape.colors[t.z], uv);
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return interpolate_quad(shape.colors[q.x], shape.colors[q.y],
        shape.colors[q.z], shape.colors[q.w], uv);
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return interpolate_line(shape.colors[l.x], shape.colors[l.y], uv.x);
  } else if (!shape.points.empty()) {
    return shape.colors[shape.points[element]];
  } else {
    return {0, 0, 0, 0};
  }
}

// Evaluate material
material_point eval_material(const scene_scene& scene,
    const scene_instance& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::glass ||
      material.type == material_type::volume ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::metallic ||
      point.type == material_type::plastic) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// check if an instance is volumetric
bool is_volumetric(const scene_scene& scene, const scene_instance& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENVIRONMENT PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate environment color.
vec3f eval_environment(const scene_scene& scene,
    const scene_environment& environment, const vec3f& direction) {
  auto wl       = transform_direction(inverse(environment.frame), direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_scene& scene, const vec3f& direction) {
  auto emission = zero3f;
  for (auto environment : scene.environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Add missing cameras.
void add_camera(scene_scene& scene) {
  scene.camera_names.emplace_back("camera");
  auto& camera        = scene.cameras.emplace_back();
  camera.orthographic = false;
  camera.film         = 0.036;
  camera.aspect       = (float)16 / (float)9;
  camera.aperture     = 0;
  camera.lens         = 0.050;
  auto bbox           = compute_bounds(scene);
  auto center         = (bbox.max + bbox.min) / 2;
  auto bbox_radius    = length(bbox.max - bbox.min) / 2;
  auto camera_dir     = vec3f{0, 0, 1};
  auto camera_dist = bbox_radius * camera.lens / (camera.film / camera.aspect);
  camera_dist *= 2.0f;  // correction for tracer camera implementation
  auto from    = camera_dir * camera_dist + center;
  auto to      = center;
  auto up      = vec3f{0, 1, 0};
  camera.frame = lookat_frame(from, to, up);
  camera.focus = length(from - to);
}

// Add a sky environment
void add_sky(scene_scene& scene, float sun_angle) {
  scene.texture_names.emplace_back("sky");
  auto& texture = scene.textures.emplace_back();
  texture       = make_sunsky(1024, 512, sun_angle);
  scene.environment_names.emplace_back("sky");
  auto& environment        = scene.environments.emplace_back();
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
}

// get named camera or default if camera is empty
camera_handle find_camera(const scene_scene& scene, const string& name) {
  if (scene.cameras.empty()) return invalid_handle;
  if (scene.camera_names.empty()) return 0;
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == name) return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "default") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera0") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera1") return idx;
  }
  return 0;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_scene& scene) {
  auto shape_bbox = vector<bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto& shape : scene.shapes) {
    auto& sbvh = shape_bbox.emplace_back();
    for (auto p : shape.positions) sbvh = merge(sbvh, p);
  }
  for (auto& instance : scene.instances) {
    auto& sbvh = shape_bbox[instance.shape];
    bbox       = merge(bbox, transform_bbox(instance.frame, sbvh));
  }
  return bbox;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

void tesselate_subdiv(
    scene_shape& shape, scene_subdiv& subdiv_, const scene_scene& scene) {
  auto subdiv = subdiv_;

  if (subdiv.subdivisions > 0) {
    if (subdiv.catmullclark) {
      std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_catmullclark(
          subdiv.quadstexcoord, subdiv.texcoords, subdiv.subdivisions, true);
      std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_catmullclark(
          subdiv.quadsnorm, subdiv.normals, subdiv.subdivisions, true);
      std::tie(subdiv.quadspos, subdiv.positions) = subdivide_catmullclark(
          subdiv.quadspos, subdiv.positions, subdiv.subdivisions);
    } else {
      std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_quads(
          subdiv.quadstexcoord, subdiv.texcoords, subdiv.subdivisions);
      std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_quads(
          subdiv.quadsnorm, subdiv.normals, subdiv.subdivisions);
      std::tie(subdiv.quadspos, subdiv.positions) = subdivide_quads(
          subdiv.quadspos, subdiv.positions, subdiv.subdivisions);
    }
    if (subdiv.smooth) {
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    } else {
      subdiv.normals   = {};
      subdiv.quadsnorm = {};
    }
  }

  if (subdiv.displacement != 0 && subdiv.displacement_tex != invalid_handle) {
    if (subdiv.texcoords.empty())
      throw std::runtime_error("missing texture coordinates");

    // facevarying case
    auto offset = vector<float>(subdiv.positions.size(), 0);
    auto count  = vector<int>(subdiv.positions.size(), 0);
    for (auto fid = 0; fid < subdiv.quadspos.size(); fid++) {
      auto qpos = subdiv.quadspos[fid];
      auto qtxt = subdiv.quadstexcoord[fid];
      for (auto i = 0; i < 4; i++) {
        auto& displacement_tex = scene.textures[subdiv.displacement_tex];
        auto  disp             = mean(
            eval_texture(displacement_tex, subdiv.texcoords[qtxt[i]], false));
        if (!displacement_tex.pixelsb.empty()) disp -= 0.5f;
        offset[qpos[i]] += subdiv.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = quads_normals(subdiv.quadspos, subdiv.positions);
    for (auto vid = 0; vid < subdiv.positions.size(); vid++) {
      subdiv.positions[vid] += normals[vid] * offset[vid] / count[vid];
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      subdiv.quadsnorm = subdiv.quadspos;
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
    }
  }

  shape = {};
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, subdiv.quadspos, subdiv.quadsnorm, subdiv.quadstexcoord,
      subdiv.positions, subdiv.normals, subdiv.texcoords);
}

void tesselate_shapes(
    scene_scene& scene, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, (int)scene.subdivs.size() + 1};
  if (progress_cb) progress_cb("tesselate subdivs", progress.x++, progress.y);

  // tesselate shapes
  for (auto& subdiv : scene.subdivs) {
    if (progress_cb) progress_cb("tesselate subdiv", progress.x++, progress.y);
    tesselate_subdiv(scene.shapes[subdiv.shape], subdiv, scene);
  }

  // done
  if (progress_cb) progress_cb("tesselate subdivs", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> scene_stats(const scene_scene& scene, bool verbose) {
  auto accumulate = [](const auto& values, const auto& func) -> size_t {
    auto sum = (size_t)0;
    for (auto& value : values) sum += func(value);
    return sum;
  };
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = compute_bounds(scene);

  auto stats = vector<string>{};
  stats.push_back("cameras:      " + format(scene.cameras.size()));
  stats.push_back("shapes:       " + format(scene.shapes.size()));
  stats.push_back("environments: " + format(scene.environments.size()));
  stats.push_back("textures:     " + format(scene.textures.size()));
  stats.push_back(
      "points:       " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.points.size(); })));
  stats.push_back(
      "lines:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.lines.size(); })));
  stats.push_back("triangles:    " +
                  format(accumulate(scene.shapes,
                      [](auto& shape) { return shape.triangles.size(); })));
  stats.push_back(
      "quads:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.quads.size(); })));
  stats.push_back("fvquads:     " +
                  format(accumulate(scene.subdivs,
                      [](auto& subdiv) { return subdiv.quadspos.size(); })));
  stats.push_back("texels4b:     " +
                  format(accumulate(scene.textures,
                      [](auto& texture) { return texture.pixelsb.size(); })));
  stats.push_back("texels4f:     " +
                  format(accumulate(scene.textures,
                      [](auto& texture) { return texture.pixelsf.size(); })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene.
vector<string> scene_validation(const scene_scene& scene, bool notextures) {
  auto errs        = vector<string>();
  auto check_names = [&errs](const vector<string>& names, const string& base) {
    auto used = unordered_map<string, int>();
    used.reserve(names.size());
    for (auto& name : names) used[name] += 1;
    for (auto& [name, used] : used) {
      if (name.empty()) {
        errs.push_back("empty " + base + " name");
      } else if (used > 1) {
        errs.push_back("duplicated " + base + " name " + name);
      }
    }
  };
  auto check_empty_textures = [&errs](const scene_scene& scene) {
    for (auto idx = 0; idx < (int)scene.textures.size(); idx++) {
      auto& texture = scene.textures[idx];
      if (texture.pixelsf.empty() && texture.pixelsb.empty()) {
        errs.push_back("empty texture " + scene.texture_names[idx]);
      }
    }
  };

  check_names(scene.camera_names, "camera");
  check_names(scene.shape_names, "shape");
  check_names(scene.material_names, "material");
  check_names(scene.instance_names, "instance");
  check_names(scene.texture_names, "texture");
  check_names(scene.environment_names, "environment");
  if (!notextures) check_empty_textures(scene);

  return errs;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(
    image_data& normalmap, const image_data& bumpmap, float scale) {
  auto width = bumpmap.width, height = bumpmap.height;
  if (normalmap.width != bumpmap.width || normalmap.height != bumpmap.height) {
    normalmap = make_image(
        width, height, bumpmap.linear, !bumpmap.pixelsf.empty());
  }
  auto dx = 1.0f / width, dy = 1.0f / height;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      auto i1 = (i + 1) % width, j1 = (j + 1) % height;
      auto p00 = get_pixel(bumpmap, i, j), p10 = get_pixel(bumpmap, i1, j),
           p01    = get_pixel(bumpmap, i, j1);
      auto g00    = (p00.x + p00.y + p00.z) / 3;
      auto g01    = (p01.x + p01.y + p01.z) / 3;
      auto g10    = (p10.x + p10.y + p10.z) / 3;
      auto normal = vec3f{
          scale * (g00 - g10) / dx, scale * (g00 - g01) / dy, 1.0f};
      normal.y = -normal.y;  // make green pointing up, even if y axis
                             // points down
      normal = normalize(normal) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
      set_pixel(normalmap, i, j, {normal.x, normal.y, normal.z, 1});
    }
  }
}
image_data bump_to_normal(const image_data& bumpmap, float scale) {
  auto normalmap = make_image(
      bumpmap.width, bumpmap.height, bumpmap.linear, !bumpmap.pixelsf.empty());
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

template <typename Shader>
static image_data make_proc_image(
    int width, int height, bool linear, bool as_byte, Shader&& shader) {
  auto image = make_image(width, height, linear, as_byte);
  auto scale = 1.0f / max(width, height);
  if (as_byte) {
    for (auto j = 0; j < height; j++) {
      for (auto i = 0; i < width; i++) {
        auto uv                      = vec2f{i * scale, j * scale};
        image.pixelsb[j * width + i] = float_to_byte(shader(uv));
      }
    }
  } else {
    for (auto j = 0; j < height; j++) {
      for (auto i = 0; i < width; i++) {
        auto uv                      = vec2f{i * scale, j * scale};
        image.pixelsf[j * width + i] = shader(uv);
      }
    }
  }
  return image;
}

// Make an image
image_data make_grid(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto thick = 0.01f / 2;
    auto c     = uv.x <= thick || uv.x >= 1 - thick || uv.y <= thick ||
             uv.y >= 1 - thick ||
             (uv.x >= 0.5f - thick && uv.x <= 0.5f + thick) ||
             (uv.y >= 0.5f - thick && uv.y <= 0.5f + thick);
    return c ? color0 : color1;
  });
}

image_data make_checker(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

image_data make_bumps(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto thick  = 0.125f;
    auto center = vec2f{
        uv.x <= 0.5f ? 0.25f : 0.75f,
        uv.y <= 0.5f ? 0.25f : 0.75f,
    };
    auto dist = clamp(length(uv - center), 0.0f, thick) / thick;
    auto val  = uv.x <= 0.5f != uv.y <= 0.5f ? (1 + sqrt(1 - dist)) / 2
                                             : (dist * dist) / 2;
    return lerp(color0, color1, val);
  });
}

image_data make_ramp(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

image_data make_gammaramp(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    if (uv.y < 1 / 3.0f) {
      return lerp(color0, color1, pow(uv.x, 2.2f));
    } else if (uv.y < 2 / 3.0f) {
      return lerp(color0, color1, uv.x);
    } else {
      return lerp(color0, color1, pow(uv.x, 1 / 2.2f));
    }
  });
}

image_data make_uvramp(int width, int height, float scale) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

image_data make_uvgrid(int width, int height, float scale, bool colored) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    uv.y     = 1 - uv.y;
    auto hsv = zero3f;
    hsv.x    = (clamp((int)(uv.x * 8), 0, 7) +
                (clamp((int)(uv.y * 8), 0, 7) + 5) % 8 * 8) /
            64.0f;
    auto vuv = uv * 4;
    vuv -= vec2f{(float)(int)vuv.x, (float)(int)vuv.y};
    auto vc  = vuv.x <= 0.5f != vuv.y <= 0.5f;
    hsv.z    = vc ? 0.5f - 0.05f : 0.5f + 0.05f;
    auto suv = uv * 16;
    suv -= vec2f{(float)(int)suv.x, (float)(int)suv.y};
    auto st = 0.01f / 2;
    auto sc = suv.x <= st || suv.x >= 1 - st || suv.y <= st || suv.y >= 1 - st;
    if (sc) {
      hsv.y = 0.2f;
      hsv.z = 0.8f;
    } else {
      hsv.y = 0.8f;
    }
    auto rgb = (colored) ? hsv_to_rgb(hsv) : vec3f{hsv.z, hsv.z, hsv.z};
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

image_data make_blackbodyramp(
    int width, int height, float scale, float from, float to) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

image_data make_colormapramp(int width, int height, float scale) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = zero3f;
    if (uv.y < 0.25) {
      rgb = colormap(uv.x, colormap_type::viridis);
    } else if (uv.y < 0.50) {
      rgb = colormap(uv.x, colormap_type::plasma);
    } else if (uv.y < 0.75) {
      rgb = colormap(uv.x, colormap_type::magma);
    } else {
      rgb = colormap(uv.x, colormap_type::inferno);
    }
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

image_data make_noisemap(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_fbmmap(int width, int height, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_turbulencemap(int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_ridgemap(int width, int height, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
image_data add_border(
    const image_data& image, float width, const vec4f& color) {
  auto result = image;
  auto scale  = 1.0f / max(image.width, image.height);
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto uv = vec2f{i * scale, j * scale};
      if (uv.x < width || uv.y < width || uv.x > image.width * scale - width ||
          uv.y > image.height * scale - width) {
        set_pixel(result, i, j, color);
      }
    }
  }
  return result;
}

// Implementation of sunsky modified heavily from pbrt
image_data make_sunsky(int width, int height, float theta_sun, float turbidity,
    bool has_sun, float sun_intensity, float sun_radius,
    const vec3f& ground_albedo) {
  auto zenith_xyY = vec3f{
      (+0.00165f * pow(theta_sun, 3.f) - 0.00374f * pow(theta_sun, 2.f) +
          0.00208f * theta_sun + 0.00000f) *
              pow(turbidity, 2.f) +
          (-0.02902f * pow(theta_sun, 3.f) + 0.06377f * pow(theta_sun, 2.f) -
              0.03202f * theta_sun + 0.00394f) *
              turbidity +
          (+0.11693f * pow(theta_sun, 3.f) - 0.21196f * pow(theta_sun, 2.f) +
              0.06052f * theta_sun + 0.25885f),
      (+0.00275f * pow(theta_sun, 3.f) - 0.00610f * pow(theta_sun, 2.f) +
          0.00316f * theta_sun + 0.00000f) *
              pow(turbidity, 2.f) +
          (-0.04214f * pow(theta_sun, 3.f) + 0.08970f * pow(theta_sun, 2.f) -
              0.04153f * theta_sun + 0.00515f) *
              turbidity +
          (+0.15346f * pow(theta_sun, 3.f) - 0.26756f * pow(theta_sun, 2.f) +
              0.06669f * theta_sun + 0.26688f),
      1000 * (4.0453f * turbidity - 4.9710f) *
              tan((4.0f / 9.0f - turbidity / 120.0f) * (pif - 2 * theta_sun)) -
          .2155f * turbidity + 2.4192f,
  };

  auto perez_A_xyY = vec3f{-0.01925f * turbidity - 0.25922f,
      -0.01669f * turbidity - 0.26078f, +0.17872f * turbidity - 1.46303f};
  auto perez_B_xyY = vec3f{-0.06651f * turbidity + 0.00081f,
      -0.09495f * turbidity + 0.00921f, -0.35540f * turbidity + 0.42749f};
  auto perez_C_xyY = vec3f{-0.00041f * turbidity + 0.21247f,
      -0.00792f * turbidity + 0.21023f, -0.02266f * turbidity + 5.32505f};
  auto perez_D_xyY = vec3f{-0.06409f * turbidity - 0.89887f,
      -0.04405f * turbidity - 1.65369f, +0.12064f * turbidity - 2.57705f};
  auto perez_E_xyY = vec3f{-0.00325f * turbidity + 0.04517f,
      -0.01092f * turbidity + 0.05291f, -0.06696f * turbidity + 0.37027f};

  auto perez_f = [](vec3f A, vec3f B, vec3f C, vec3f D, vec3f E, float theta,
                     float gamma, float theta_sun, vec3f zenith) -> vec3f {
    auto num = ((1 + A * exp(B / cos(theta))) *
                (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
    auto den = ((1 + A * exp(B)) * (1 + C * exp(D * theta_sun) +
                                       E * cos(theta_sun) * cos(theta_sun)));
    return zenith * num / den;
  };

  auto sky = [&perez_f, perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                 perez_E_xyY, zenith_xyY](
                 float theta, float gamma, float theta_sun) -> vec3f {
    return xyz_to_rgb(xyY_to_xyz(
               perez_f(perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                   perez_E_xyY, theta, gamma, theta_sun, zenith_xyY))) /
           10000;
  };

  // compute sun luminance
  auto sun_ko     = vec3f{0.48f, 0.75f, 0.14f};
  auto sun_kg     = vec3f{0.1f, 0.0f, 0.0f};
  auto sun_kwa    = vec3f{0.02f, 0.0f, 0.0f};
  auto sun_sol    = vec3f{20000.0f, 27000.0f, 30000.0f};
  auto sun_lambda = vec3f{680, 530, 480};
  auto sun_beta   = 0.04608365822050f * turbidity - 0.04586025928522f;
  auto sun_m      = 1.0f /
               (cos(theta_sun) + 0.000940f * pow(1.6386f - theta_sun, -1.253f));

  auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda / 1000, -4.08f));
  auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda / 1000, -1.3f));
  auto tauO = exp(-sun_m * sun_ko * .35f);
  auto tauG = exp(
      -1.41f * sun_kg * sun_m / pow(1 + 118.93f * sun_kg * sun_m, 0.45f));
  auto tauWA  = exp(-0.2385f * sun_kwa * 2.0f * sun_m /
                   pow(1 + 20.07f * sun_kwa * 2.0f * sun_m, 0.45f));
  auto sun_le = sun_sol * tauR * tauA * tauO * tauG * tauWA * 10000;

  // rescale by user
  sun_le *= sun_intensity;

  // sun scale from Wikipedia scaled by user quantity and rescaled to at
  // the minimum 5 pixel diamater
  auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
  sun_angular_radius *= sun_radius;
  sun_angular_radius = max(sun_angular_radius, 2 * pif / height);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000 : zero3f;
  };

  // Make the sun sky image
  auto img          = make_image(width, height, true, false);
  auto sky_integral = 0.0f, sun_integral = 0.0f;
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < width; i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / width);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      sky_integral += mean(sky_col) * sin(theta);
      sun_integral += mean(sun_col) * sin(theta);
      auto col                   = sky_col + sun_col;
      img.pixelsf[j * width + i] = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != zero3f) {
    auto ground = zero3f;
    for (auto j = 0; j < height / 2; j++) {
      auto theta = pif * ((j + 0.5f) / height);
      for (int i = 0; i < width; i++) {
        auto pxl   = img.pixelsf[j * width + i];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (width * height);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        img.pixelsf[j * width + i] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        img.pixelsf[j * width + i] = {0, 0, 0, 1};
      }
    }
  }

  // done
  return img;
}

// Make an image of multiple lights.
image_data make_lights(int width, int height, const vec3f& le, int nlights,
    float langle, float lwidth, float lheight) {
  auto img = make_image(width, height, true, false);
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < width; i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / width);
      auto inlight = false;
      for (auto l = 0; l < nlights; l++) {
        auto lphi = 2 * pif * (l + 0.5f) / nlights;
        inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
      }
      img.pixelsf[j * width + i] = {le.x, le.y, le.z, 1};
    }
  }
  return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a plane.
shape_data make_rect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_rect(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
shape_data make_bulged_rect(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = shape_data{};
  make_bulged_rect(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a plane in the xz plane.
shape_data make_recty(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_recty(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
shape_data make_bulged_recty(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = shape_data{};
  make_bulged_recty(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Make a box.
shape_data make_box(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = shape_data{};
  make_box(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
shape_data make_rounded_box(const vec3i& steps, const vec3f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = shape_data{};
  make_rounded_box(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a quad stack
shape_data make_rect_stack(
    const vec3i& steps, const vec3f& scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_rect_stack(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a floor.
shape_data make_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
shape_data make_bent_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale, float bent) {
  auto shape = shape_data{};
  make_bent_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, bent);
  return shape;
}

// Make a sphere.
shape_data make_sphere(int steps, float scale, float uvscale) {
  auto shape = shape_data{};
  make_sphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere.
shape_data make_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_uvsphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere with slipped caps.
shape_data make_capped_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale, float height) {
  auto shape = shape_data{};
  make_capped_uvsphere(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, height);
  return shape;
}
// Make a disk
shape_data make_disk(int steps, float scale, float uvscale) {
  auto shape = shape_data{};
  make_disk(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}

// Make a bulged disk
shape_data make_bulged_disk(
    int steps, float scale, float uvscale, float height) {
  auto shape = shape_data{};
  make_bulged_disk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, height);
  return shape;
}

// Make a uv disk
shape_data make_uvdisk(const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = shape_data{};
  make_uvdisk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a uv cylinder
shape_data make_uvcylinder(
    const vec3i& steps, const vec2f& scale, const vec3f& uvscale) {
  auto shape = shape_data{};
  make_uvcylinder(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(const vec3i& steps, const vec2f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = shape_data{};
  make_rounded_uvcylinder(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
shape_data make_lines(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, const vec2f& rad) {
  auto shape = shape_data{};
  make_lines(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, steps, scale, uvscale, rad);
  return shape;
}

// Make point primitives. Returns points, pos, norm, texcoord, radius.
shape_data make_point(float radius) {
  auto shape = shape_data{};
  make_point(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, radius);
  return shape;
}

shape_data make_points(int num, float uvscale, float radius) {
  auto shape = shape_data{};
  make_points(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, num, uvscale, radius);
  return shape;
}

shape_data make_random_points(
    int num, const vec3f& size, float uvscale, float radius, uint64_t seed) {
  auto shape = shape_data{};
  make_random_points(shape.points, shape.positions, shape.normals,
      shape.texcoords, shape.radius, num, size, uvscale, radius, seed);
  return shape;
}

// Make a facevarying rect
fvshape_data make_fvrect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = fvshape_data{};
  make_fvrect(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying box
fvshape_data make_fvbox(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = fvshape_data{};
  make_fvbox(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying sphere
fvshape_data make_fvsphere(int steps, float scale, float uvscale) {
  auto shape = fvshape_data{};
  make_fvsphere(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Predefined meshes
shape_data make_monkey(float scale) {
  auto shape = shape_data{};
  make_monkey(shape.quads, shape.positions, scale);
  return shape;
}
shape_data make_quad(float scale) {
  auto shape = shape_data{};
  make_quad(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_quady(float scale) {
  auto shape = shape_data{};
  make_quady(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_cube(float scale) {
  auto shape = shape_data{};
  make_cube(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
fvshape_data make_fvcube(float scale) {
  auto shape = fvshape_data{};
  make_fvcube(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_geosphere(float scale) {
  auto shape = shape_data{};
  make_geosphere(shape.triangles, shape.positions, scale);
  return shape;
}

// Make a hair ball around a shape
shape_data make_hair(const shape_data& base, const vec2i& steps,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
  auto shape = shape_data{};
  make_hair(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, base.triangles, base.quads, base.positions, base.normals,
      base.texcoords, steps, len, rad, noise, clump, rotation, seed);
  return shape;
}

// Make a heightfield mesh.
shape_data make_heightfield(const vec2i& size, const vector<float>& height) {
  auto shape = shape_data{};
  make_heightfield(shape.quads, shape.positions, shape.normals, shape.texcoords,
      size, height);
  return shape;
}
shape_data make_heightfield(const vec2i& size, const vector<vec4f>& color) {
  auto shape = shape_data{};
  make_heightfield(shape.quads, shape.positions, shape.normals, shape.texcoords,
      size, color);
  return shape;
}

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res and without texcoords and normals.
shape_data points_to_spheres(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  points_to_spheres(shape.quads, shape.positions, shape.normals,
      shape.texcoords, vertices, steps, scale);
  return shape;
}
shape_data lines_to_cylinders(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  lines_to_cylinders(shape.quads, shape.positions, shape.normals,
      shape.texcoords, vertices, steps, scale);
  return shape;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox(scene_scene& scene) {
  scene.asset.name = "cornellbox";

  auto& camera    = scene.cameras.emplace_back();
  camera.frame    = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 3.9}};
  camera.lens     = 0.035;
  camera.aperture = 0.0;
  camera.focus    = 3.9;
  camera.film     = 0.024;
  camera.aspect   = 1;

  auto& floor_shape       = scene.shapes.emplace_back();
  floor_shape.positions   = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& floor_material    = scene.materials.emplace_back();
  floor_material.color    = {0.725, 0.71, 0.68};
  auto& floor_instance    = scene.instances.emplace_back();
  floor_instance.shape    = (int)scene.shapes.size() - 1;
  floor_instance.material = (int)scene.materials.size() - 1;

  auto& ceiling_shape       = scene.shapes.emplace_back();
  ceiling_shape.positions   = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& ceiling_material    = scene.materials.emplace_back();
  ceiling_material.color    = {0.725, 0.71, 0.68};
  auto& ceiling_instance    = scene.instances.emplace_back();
  ceiling_instance.shape    = (int)scene.shapes.size() - 1;
  ceiling_instance.material = (int)scene.materials.size() - 1;

  auto& backwall_shape     = scene.shapes.emplace_back();
  backwall_shape.positions = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& backwall_material  = scene.materials.emplace_back();
  backwall_material.color  = {0.725, 0.71, 0.68};
  auto& backwall_instance  = scene.instances.emplace_back();
  backwall_instance.shape  = (int)scene.shapes.size() - 1;
  backwall_instance.material = (int)scene.materials.size() - 1;

  auto& rightwall_shape       = scene.shapes.emplace_back();
  rightwall_shape.positions   = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& rightwall_material    = scene.materials.emplace_back();
  rightwall_material.color    = {0.14, 0.45, 0.091};
  auto& rightwall_instance    = scene.instances.emplace_back();
  rightwall_instance.shape    = (int)scene.shapes.size() - 1;
  rightwall_instance.material = (int)scene.materials.size() - 1;

  auto& leftwall_shape     = scene.shapes.emplace_back();
  leftwall_shape.positions = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& leftwall_material  = scene.materials.emplace_back();
  leftwall_material.color  = {0.63, 0.065, 0.05};
  auto& leftwall_instance  = scene.instances.emplace_back();
  leftwall_instance.shape  = (int)scene.shapes.size() - 1;
  leftwall_instance.material = (int)scene.materials.size() - 1;

  auto& shortbox_shape       = scene.shapes.emplace_back();
  shortbox_shape.positions   = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  shortbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& shortbox_material    = scene.materials.emplace_back();
  shortbox_material.color    = {0.725, 0.71, 0.68};
  auto& shortbox_instance    = scene.instances.emplace_back();
  shortbox_instance.shape    = (int)scene.shapes.size() - 1;
  shortbox_instance.material = (int)scene.materials.size() - 1;

  auto& tallbox_shape       = scene.shapes.emplace_back();
  tallbox_shape.positions   = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
      {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
      {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
      {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
      {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
      {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
      {-0.71, 0.0, -0.49}};
  tallbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& tallbox_material    = scene.materials.emplace_back();
  tallbox_material.color    = {0.725, 0.71, 0.68};
  auto& tallbox_instance    = scene.instances.emplace_back();
  tallbox_instance.shape    = (int)scene.shapes.size() - 1;
  tallbox_instance.material = (int)scene.materials.size() - 1;

  auto& light_shape       = scene.shapes.emplace_back();
  light_shape.positions   = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
      {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  light_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& light_material    = scene.materials.emplace_back();
  light_material.emission = {17, 12, 4};
  auto& light_instance    = scene.instances.emplace_back();
  light_instance.shape    = (int)scene.shapes.size() - 1;
  light_instance.material = (int)scene.materials.size() - 1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Find the first keyframe value that is greater than the argument.
inline int keyframe_index(const vector<float>& times, const float& time) {
  for (auto i = 0; i < times.size(); i++)
    if (times[i] > time) return i;
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
