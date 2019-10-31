//
// # Yocto/Image: Tiny imaging Library mostly for rendering and color support
//
//
// Yocto/Image is a collection of image utilities useful when writing rendering
// algorithms. These include a simple image data structure, color conversion
// utilities and tone mapping. We provinde loading and saving functionality for
// images and support PNG, JPG, TGA, BMP, HDR, EXR formats.
//
// This library depends on stb_image.h, stb_image_write.h, stb_image_resize.h,
// tinyexr.h for the IO features. If thoese are not needed, it can be safely
// used without dependencies.
//
//
// ## Image Utilities
//
// Yocto/Image supports a very small set of color and image utilities including
// color utilities, example image creation, tone mapping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// global illumination renderer, rather than the need of generic image editing.
// We support 4-channels float images (assumed to be in linear color) and
// 4-channels byte images (assumed to be in sRGB).
//
//
// 1. store images using the image<T> structure
// 2. load and save images with `load_image()` and `save_image()`
// 3. resize images with `resize()`
// 4. tonemap images with `tonemap()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_proc_image()` functions
// 6. create procedural sun-sky images with `make_sunsky()`
// 7. many color conversion functions are available in the code below
//
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
//
//  LICENSE for blackbody code
//
// Copyright (c) 2015 Neil Bartlett
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
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//

#ifndef _YOCTO_IMAGE_H_
#define _YOCTO_IMAGE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates a color image at a point `uv`.
vec4f eval_image(const image<vec4f>& img, const vec2f& uv,
    bool no_interpolation, bool clamp_to_edge);
vec4f eval_image(const image<vec4b>& img, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt);
image<vec4b> float_to_byte(const image<vec4f>& fl);
void         byte_to_float(image<vec4f>& fl, const image<vec4b>& bt);
void         float_to_byte(image<vec4b>& bt, const image<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_rgb(const image<vec4f>& srgb);
image<vec4f> rgb_to_srgb(const image<vec4f>& rgb);
image<vec4f> srgb_to_rgb(const image<vec4b>& srgb);
image<vec4b> rgb_to_srgbb(const image<vec4f>& rgb);
void         srgb_to_rgb(image<vec4f>& rgb, const image<vec4f>& srgb);
void         rgb_to_srgb(image<vec4f>& srgb, const image<vec4f>& rgb);

// Tone mapping params
struct tonemap_params {
  float exposure    = 0;
  vec3f tint        = {1, 1, 1};
  float contrast    = 0.5;
  float logcontrast = 0.5;
  float saturation  = 0.5;
  bool  filmic      = false;
  bool  srgb        = true;
};

// Apply tone mapping
vec3f tonemap(const vec3f& hdr, const tonemap_params& params);
vec4f tonemap(const vec4f& hdr, const tonemap_params& params);

// Apply tone mapping
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, const tonemap_params& params);
image<vec4b> tonemap_imageb(
    const image<vec4f>& hdr, const tonemap_params& params);

// minimal color grading
struct colorgrade_params {
  float contrast         = 0.5;
  float shadows          = 0.5;
  float midtones         = 0.5;
  float highlights       = 0.5;
  vec3f shadows_color    = {1, 1, 1};
  vec3f midtones_color   = {1, 1, 1};
  vec3f highlights_color = {1, 1, 1};
};

vec3f colorgrade(const vec3f& ldr, const colorgrade_params& params);
vec4f colorgrade(const vec4f& ldr, const colorgrade_params& params);

// color grade an image region
image<vec4f> colorgrade_image(
    const image<vec4f>& img, const colorgrade_params& params);
void colorgrade_region(image<vec4f>& corrected, const image<vec4f>& img,
    const image_region& region, const colorgrade_params& params);

// determine white balance colors
vec3f compute_white_balance(const image<vec4f>& img);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size);
image<vec4b> resize_image(const image<vec4b>& img, const vec2i& size);
void         resize_image(
            image<vec4f>& res, const image<vec4f>& img, const vec2i& size);
void resize_image(
    image<vec4b>& res, const image<vec4b>& img, const vec2i& size);

// Compute the difference between two images
image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool disply_diff);
void image_difference(image<vec4f>& diff, const image<vec4f>& a,
    const image<vec4f>& b, bool disply_diff);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Result of io operations
struct imageio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
image<vec4f>   load_image(const string& filename);
imageio_status load_image(const string& filename, image<vec4f>& img);
imageio_status save_image(const string& filename, const image<vec4f>& img);
image<vec4b>   load_imageb(const string& filename);
imageio_status load_imageb(const string& filename, image<vec4b>& img);
imageio_status save_imageb(const string& filename, const image<vec4b>& img);

// Loads/saves a 4 channels float/byte image tonemapped
imageio_status save_image_tonemapped(const string& filename, const image<vec4f>& img, const tonemap_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Parameters for make_proc_image
struct proc_image_params {
  // clang-format off
  enum struct type_t {
    grid, checker, bumps, ramp, gammaramp, uvramp, uvgrid, blackbody, noise,
    turbulence, fbm, ridge };
  // clang-format on
  type_t type    = type_t::grid;
  vec2i  size    = {1024, 1024};
  float  scale   = 1;
  vec4f  color0  = {0, 0, 0, 1};
  vec4f  color1  = {1, 1, 1, 1};
  vec4f  noise   = {2, 0.5, 8, 1};  // lacunarity, gain, octaves, offset
  float  borderw = 0;
  vec4f  borderc = {0, 0, 0, 1};
};

// Make an image
image<vec4f> make_proc_image(const proc_image_params& params);
void make_proc_image(image<vec4f>& img, const proc_image_params& params);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
image<vec4f> make_sunsky(const vec2i& size, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
    float sun_radius = 1, const vec3f& ground_albedo = {0.2, 0.2, 0.2});
void         make_sunsky(image<vec4f>& img, const vec2i& size, float sun_angle,
            float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
            float sun_radius = 1, const vec3f& ground_albedo = {0.2, 0.2, 0.2});
// Make an image of multiple lights.
image<vec4f> make_lights(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);
void         make_lights(image<vec4f>& img, const vec2i& size,
            const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
            float lwidth = pif / 16, float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
image<vec4f> bump_to_normal(const image<vec4f>& img, float scale = 1);
void         bump_to_normal(
            image<vec4f>& norm, const image<vec4f>& img, float scale = 1);

// Add a border to an image
image<vec4f> add_border(const image<vec4f>& img, int width, const vec4f& color);
void         add_border(
            image<vec4f>& bordered, image<vec4f>& img, int width, const vec4f& color);

// Make logo images. Image is resized to proper size.
image<vec4b> make_logo(const string& name);
void         make_logo(image<vec4f>& img, const string& name);
void         make_logo(image<vec4b>& img, const string& name);
image<vec4f> add_logo(
    const image<vec4f>& img, const string& name = "logo-medium");
image<vec4b> add_logo(
    const image<vec4b>& img, const string& name = "logo-medium");
void add_logo(image<vec4f>& with_logo, const image<vec4f>& img,
    const string& name = "logo-medium");
void add_logo(image<vec4b>& with_logo, const image<vec4b>& img,
    const string& name = "logo-medium");

// Make an image preset, useful for testing. See implementation for types.
image<vec4f> make_image_preset(const string& type);
bool         make_image_preset(image<vec4f>& img, const string& type);
bool         make_image_preset(image<vec4b>& img, const string& type);
bool         make_image_preset(
            image<vec4f>& hdr, image<vec4b>& ldr, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates a color image at a point `uv`.
float eval_volume(const image<float>& img, const vec3f& uvw,
    bool no_interpolation = false, bool clamp_to_edge = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1 channel volume.
void load_volume(const string& filename, volume<float>& vol);
void save_volume(const string& filename, const volume<float>& vol);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE VOLUMES
// -----------------------------------------------------------------------------
namespace yocto {

// make a simple example volume
void make_voltest(volume<float>& vol, const vec3i& size, float scale = 10,
    float exponent = 6);
void make_volume_preset(volume<float>& vol, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// RGB color spaces
enum struct color_space {
  rgb,         // default linear space (srgb linear)
  srgb,        // srgb color space (non-linear)
  adobe,       // Adobe rgb color space (non-linear)
  prophoto,    // ProPhoto Kodak rgb color space (non-linear)
  rec709,      // hdtv color space (non-linear)
  rec2020,     // uhtv color space (non-linear)
  rec2100pq,   // hdr color space with perceptual quantizer (non-linear)
  rec2100hlg,  // hdr color space with hybrid log gamma (non-linear)
  aces2065,    // ACES storage format (linear)
  acescg,      // ACES CG computation (linear)
  acescc,      // ACES color correction (non-linear)
  acescct,     // ACES color correction 2 (non-linear)
  p3dci,       // P3 DCI (non-linear)
  p3d60,       // P3 variation for D60 (non-linear)
  p3d65,       // P3 variation for D65 (non-linear)
  p3display,   // Apple display P3
};

// Conversion between rgb color spaces
vec3f color_to_xyz(const vec3f& col, color_space from);
vec3f xyz_to_color(const vec3f& xyz, color_space to);

}  // namespace yocto

#endif
