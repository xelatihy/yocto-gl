//
// # Yocto/Pbrt: Tiny library for Pbrt parsing
//
// Yocto/Pbrt is a simple pbrt parser that works with callbacks.
// We make no attempt to provide a simple interface for pbrt but just the
// low level parsing code.
//
// Error reporting is done through exceptions using the `io_error`
// exception.
//
// ## Parse an pbrt file
//
// 1. define callbacks in `pbrt_callback` structure using lambda with capture
//    if desired
// 2. run the parse with `load_pbrt()`
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

#ifndef _YOCTO_PBRT_H_
#define _YOCTO_PBRT_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_utils.h"

#include <string_view>
#include <unordered_map>
#include <variant>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::function;
using std::get;
using std::holds_alternative;
using std::string_view;
using std::unordered_map;
using std::variant;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// pbrt pbrt_spectrum as rgb color
struct pbrt_spectrum3f {
    float x, y, z;

    constexpr pbrt_spectrum3f() : x{0}, y{0}, z{0} {}
    constexpr pbrt_spectrum3f(float x, float y, float z) : x{x}, y{y}, z{z} {}
    constexpr explicit pbrt_spectrum3f(float v) : x{v}, y{v}, z{v} {}
    constexpr explicit operator vec3f() const { return {x, y, z}; };

    constexpr float&       operator[](int i) { return (&x)[i]; }
    constexpr const float& operator[](int i) const { return (&x)[i]; }
};

// pbrt cameras
struct pbrt_perspective_camera {
    float  fov              = 90;
    float  frameaspectratio = -1;  // or computed from film
    float  lensradius       = 0;
    float  focaldistance    = 1e30;
    bbox2f screenwindow     = {{-1, -1}, {1, 1}};
    float  shutteropen      = 0;
    float  shutterclose     = 1;
};
struct pbrt_orthographic_camera {
    float  frameaspectratio = -1;  // or computed from film
    float  lensradius       = 0;
    float  focaldistance    = 1e30;
    bbox2f screenwindow     = {{-1, -1}, {1, 1}};
    float  shutteropen      = 0;
    float  shutterclose     = 1;
};
struct pbrt_environment_camera {
    float shutteropen  = 0;
    float shutterclose = 1;
};
struct pbrt_realistic_camera {
    string lensfile           = "";
    float  aperturediameter   = 1;
    float  focusdistance      = 10;
    bool   simpleweighting    = true;
    float  shutteropen        = 0;
    float  shutterclose       = 1;
    float  approx_focallength = 0;
};
enum struct pbrt_camera_type { perspective, orthographic, environment, realistic };
struct pbrt_camera {
    pbrt_camera_type type = pbrt_camera_type::perspective;
    pbrt_perspective_camera perspective = {};
    pbrt_orthographic_camera orthographic = {};
    pbrt_environment_camera environment = {};
    pbrt_realistic_camera realistic = {};
};

// pbrt samplers
struct pbrt_random_sampler {
    int pixelsamples = 16;
};
struct pbrt_halton_sampler {
    int pixelsamples = 16;
};
struct pbrt_sobol_sampler {
    int pixelsamples = 16;
};
struct pbrt_zerotwosequence_sampler {
    int pixelsamples = 16;
};
struct pbrt_maxmindist_sampler {
    int pixelsamples = 16;
};
struct pbrt_stratified_sampler {
    bool jitter   = true;
    int  xsamples = 2;
    int  ysamples = 2;
};
enum struct pbrt_sampler_type { random, halton, sobol, zerotwosequence, maxmindist, stratified };
struct pbrt_sampler {
    pbrt_sampler_type type = pbrt_sampler_type::random;
    pbrt_random_sampler random = {};
    pbrt_halton_sampler halton = {};
    pbrt_sobol_sampler sobol = {};
    pbrt_zerotwosequence_sampler zerotwosequence = {};
    pbrt_maxmindist_sampler maxmindist = {};
    pbrt_stratified_sampler stratified = {};

};

// pbrt film
struct pbrt_image_film {
    int    xresolution        = 640;
    int    yresolution        = 480;
    bbox2f cropwindow         = {{0, 0}, {1, 1}};
    float  scale              = 1;
    float  maxsampleluminance = float_max;
    float  diagonal           = 35;
    string filename           = "pbrt.exr";
};
enum struct pbrt_film_type { image };
struct pbrt_film {
    pbrt_film_type type = pbrt_film_type::image;
    pbrt_image_film image = {};
};

// pbrt filters
struct pbrt_box_filter {
    float xwidth = 0.5f;
    float ywidth = 0.5f;
};
struct pbrt_gaussian_filter {
    float xwidth = 2;
    float ywidth = 2;
    float alpha  = 2;
};
struct pbrt_mitchell_filter {
    float xwidth = 2;
    float ywidth = 2;
    float B      = 1.0f / 3.0f;
    float C      = 1.0f / 3.0f;
};
struct pbrt_sinc_filter {
    float xwidth = 4;
    float ywidth = 4;
    float tau    = 3;
};
struct pbrt_triangle_filter {
    float xwidth = 2;
    float ywidth = 2;
};
enum struct pbrt_filter_type { box, gaussian, mitchell, sinc, triangle };
struct pbrt_filter {
    pbrt_filter_type type = pbrt_filter_type::box;
    pbrt_box_filter box = {};
    pbrt_gaussian_filter gaussian = {};
    pbrt_mitchell_filter mitchell = {};
    pbrt_sinc_filter sinc = {};
    pbrt_triangle_filter triangle = {};
};

// pbrt integrators
struct pbrt_path_integrator {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    float  rrthreshold = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
};
struct pbrt_volpath_integrator {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    float  rrthreshold = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
};
struct pbrt_bdpt_integrator {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::power;
    bool                  visualizestrategies = false;
    bool                  visualizeweights    = false;
};
struct pbrt_directlighting_integrator {
    enum struct strategy_t { all, one };
    strategy_t strategy    = strategy_t::all;
    int        maxdepth    = 5;
    bbox2i     pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
};
struct pbrt_mlt_integrator {
    int    maxdepth             = 5;
    bbox2i pixelbounds          = {{0, 0}, {type_max<int>, type_max<int>}};
    int    bootstrapsamples     = 100000;
    int    chains               = 1000;
    int    mutationsperpixel    = 100;
    float  largestepprobability = 0.3;
    float  sigma                = 0.01;
};
struct pbrt_sppm_integrator {
    int    maxdepth            = 5;
    bbox2i pixelbounds         = {{0, 0}, {type_max<int>, type_max<int>}};
    int    iterations          = 64;
    int    photonsperiteration = -1;
    int    imagewritefrequency = pow2(31);
    float  radius              = 5;
};
struct pbrt_whitted_integrator {
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
};
enum struct pbrt_integrator_type { path, volpath, bdpt, directlighting, mlt, sppm, whitted };
struct pbrt_integrator {
    pbrt_integrator_type type = pbrt_integrator_type::path;
    pbrt_path_integrator path = {};
    pbrt_volpath_integrator volpath = {};
    pbrt_bdpt_integrator bdpt = {};
    pbrt_directlighting_integrator directlighting = {};
    pbrt_mlt_integrator mlt = {};
    pbrt_sppm_integrator sppm = {};
    pbrt_whitted_integrator whitted = {};
};

// pbrt accellerators
struct pbrt_bvh_accelerator {
    enum struct splitmethod_t { sah, equal, middle, hlbvh };
    int           maxnodeprims = 4;
    splitmethod_t splitmethod  = splitmethod_t::sah;
};
struct pbrt_kdtree_accelerator {
    int   intersectcost = 80;
    int   traversalcost = 1;
    float emptybonus    = 0.2;
    int   maxprims      = 1;
    int   maxdepth      = -1;
};
enum struct pbrt_accelerator_type { bvh, kdtree };
struct pbrt_accelerator {
    pbrt_accelerator_type type = pbrt_accelerator_type::bvh;
    pbrt_bvh_accelerator bvh = {};
    pbrt_kdtree_accelerator kdtree = {};
};

// pbrt texture or value
template <typename T>
struct pbrt_textured;

template <>
struct pbrt_textured<float> {
    float  value   = 0;
    string texture = "";
    pbrt_textured() : value{0}, texture{} {}
    pbrt_textured(float v) : value{v}, texture{} {}
};
template <>
struct pbrt_textured<pbrt_spectrum3f> {
    pbrt_spectrum3f value   = {0, 0, 0};
    string          texture = "";
    pbrt_textured() : value{0, 0, 0}, texture{} {}
    pbrt_textured(float x, float y, float z) : value{x, y, z}, texture{} {}
};

// pbrt textures
struct pbrt_constant_texture {
    pbrt_textured<pbrt_spectrum3f> value = {1, 1, 1};
};
struct pbrt_bilerp_texture {
    pbrt_textured<pbrt_spectrum3f> v00 = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> v01 = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> v10 = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> v11 = {1, 1, 1};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_checkerboard_texture {
    enum struct aamode_type { closedform, none };
    int                            dimension = 2;
    pbrt_textured<pbrt_spectrum3f> tex1      = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> tex2      = {0, 0, 0};
    aamode_type                    aamode    = aamode_type::closedform;
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_dots_texture {
    pbrt_textured<pbrt_spectrum3f> inside  = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> outside = {0, 0, 0};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_fbm_texture {
    int   octaves   = 8;
    float roughness = 0.5;
};
struct pbrt_imagemap_texture {
    enum wrap_type { repeat, black, clamp };
    string    filename      = "";
    wrap_type wrap          = wrap_type::repeat;
    float     maxanisotropy = 8;
    bool      trilinear     = false;
    float     scale         = 1;
    bool      gamma         = true;
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_marble_texture {
    int   octaves   = 8;
    float roughness = 0.5f;
    float scale     = 1;
    float variation = 0.2f;
};
struct pbrt_mix_texture {
    pbrt_textured<pbrt_spectrum3f> tex1   = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> tex2   = {1, 1, 1};
    pbrt_textured<float>           amount = 0.5f;
};
struct pbrt_scale_texture {
    pbrt_textured<pbrt_spectrum3f> tex1 = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> tex2 = {1, 1, 1};
};
struct pbrt_uv_texture {
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_windy_texture {
    // TODO: missing parameters
};
struct pbrt_wrinkled_texture {
    int   octaves   = 8;
    float roughness = 0.5;
};
using pbrt_texture = variant<pbrt_constant_texture, pbrt_bilerp_texture,
    pbrt_checkerboard_texture, pbrt_dots_texture, pbrt_fbm_texture,
    pbrt_imagemap_texture, pbrt_marble_texture, pbrt_mix_texture,
    pbrt_scale_texture, pbrt_uv_texture, pbrt_windy_texture,
    pbrt_wrinkled_texture>;

// pbrt materials
struct pbrt_matte_material {
    pbrt_textured<pbrt_spectrum3f> Kd      = {0.5, 0.5, 0.5};
    pbrt_textured<float>           sigma   = 0;
    pbrt_textured<float>           bumpmap = 0;
};
struct pbrt_mirror_material {
    pbrt_textured<pbrt_spectrum3f> Kr      = {0.9, 0.9, 0.9};
    pbrt_textured<float>           bumpmap = 0;
};
struct pbrt_plastic_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.25, 0.25, 0.25};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0.25, 0.25, 0.25};
    pbrt_textured<float>           uroughness     = 0.1;
    pbrt_textured<float>           vroughness     = 0.1;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_metal_material {
    pbrt_textured<pbrt_spectrum3f> eta = {
        0.2004376970f, 0.9240334304f, 1.1022119527f};
    pbrt_textured<pbrt_spectrum3f> k = {
        3.9129485033f, 2.4528477015f, 2.1421879552f};
    pbrt_textured<float> uroughness     = 0.01;
    pbrt_textured<float> vroughness     = 0.01;
    bool                 remaproughness = true;
    pbrt_textured<float> bumpmap        = 0;
};
struct pbrt_glass_material {
    pbrt_textured<pbrt_spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>           eta            = 1.5;
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_translucent_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.25, 0.25, 0.25};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0.25, 0.25, 0.25};
    pbrt_textured<pbrt_spectrum3f> reflect        = {0.5, 0.5, 0.5};
    pbrt_textured<pbrt_spectrum3f> transmit       = {0.5, 0.5, 0.5};
    pbrt_textured<float>           uroughness     = 0.1;
    pbrt_textured<float>           vroughness     = 0.1;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_uber_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.25, 0.25, 0.25};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0.25, 0.25, 0.25};
    pbrt_textured<pbrt_spectrum3f> Kr             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Kt             = {0, 0, 0};
    pbrt_textured<float>           uroughness     = 0.1;
    pbrt_textured<float>           vroughness     = 0.1;
    pbrt_textured<float>           eta            = 1.5;
    pbrt_textured<pbrt_spectrum3f> opacity        = {1, 1, 1};
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_disney_material {
    pbrt_textured<pbrt_spectrum3f> color           = {0.5f, 0.5f, 0.5f};
    pbrt_textured<float>           anisotropic     = 0;
    pbrt_textured<float>           clearcoat       = 0;
    pbrt_textured<float>           clearcoatgloss  = 1;
    pbrt_textured<float>           eta             = 1.5;
    pbrt_textured<float>           metallic        = 0;
    pbrt_textured<float>           uroughness      = 0.5;
    pbrt_textured<float>           vroughness      = 0.5;
    pbrt_textured<pbrt_spectrum3f> scatterdistance = {0, 0, 0};
    pbrt_textured<float>           sheen           = 0;
    pbrt_textured<float>           sheentint       = 0.5;
    pbrt_textured<float>           spectrans       = 0;
    pbrt_textured<float>           speculartint    = 0;
    bool                           thin            = false;
    pbrt_textured<pbrt_spectrum3f> difftrans       = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> flatness        = {0, 0, 0};
    bool                           remaproughness  = true;
    pbrt_textured<float>           bumpmap         = 0;
};
struct pbrt_fourier_material {
    string               bsdffile = "";
    pbrt_textured<float> bumpmap  = 0;
    variant<pbrt_plastic_material, pbrt_metal_material, pbrt_glass_material>
        approx = {};
};
struct pbrt_hair_material {
    pbrt_textured<pbrt_spectrum3f> color = {0, 0, 0};  // TODO: missing default
    pbrt_textured<pbrt_spectrum3f> sigma_a = {
        0, 0, 0};                          // TODO: missing default
    pbrt_textured<float> eumelanin   = 0;  // TODO: missing default
    pbrt_textured<float> pheomelanin = 0;  // TODO: missing default
    pbrt_textured<float> eta         = 1.55f;
    pbrt_textured<float> beta_m      = 0.3f;
    pbrt_textured<float> beta_n      = 0.3f;
    pbrt_textured<float> alpha       = 2;
    pbrt_textured<float> bumpmap     = 0;
};
struct pbrt_kdsubsurface_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.5, 0.5, 0.5};
    pbrt_textured<pbrt_spectrum3f> mfp            = {1, 1, 1};
    pbrt_textured<float>           eta            = 1.3;
    pbrt_textured<pbrt_spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_mix_material {
    pbrt_textured<pbrt_spectrum3f> amount         = {0, 0, 0};
    string                         namedmaterial1 = "";
    string                         namedmaterial2 = "";
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_substrate_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.5, 0.5, 0.5};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0.5, 0.5, 0.5};
    pbrt_textured<float>           uroughness     = 0.1;
    pbrt_textured<float>           vroughness     = 0.1;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_subsurface_material {
    string                         name           = "";
    pbrt_textured<pbrt_spectrum3f> sigma_a        = {.0011, .0024, .014};
    pbrt_textured<pbrt_spectrum3f> sigma_prime_s  = {2.55, 3.12, 3.77};
    float                          scale          = 1;
    pbrt_textured<float>           eta            = 1;
    pbrt_textured<pbrt_spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<pbrt_spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
using pbrt_material = variant<pbrt_matte_material, pbrt_mirror_material,
    pbrt_plastic_material, pbrt_metal_material, pbrt_glass_material,
    pbrt_translucent_material, pbrt_uber_material, pbrt_disney_material,
    pbrt_fourier_material, pbrt_hair_material, pbrt_kdsubsurface_material,
    pbrt_mix_material, pbrt_substrate_material, pbrt_subsurface_material>;

// pbrt shapes
struct pbrt_trianglemesh_shape {
    vector<vec3i>        indices     = {};
    vector<vec3f>        P           = {};
    vector<vec3f>        N           = {};
    vector<vec3f>        S           = {};
    vector<vec2f>        uv          = {};
    pbrt_textured<float> alpha       = 1;
    pbrt_textured<float> shadowalpha = 1;
};
struct pbrt_plymesh_shape {
    string               filename    = {};
    pbrt_textured<float> alpha       = 1;
    pbrt_textured<float> shadowalpha = 1;
};
struct pbrt_curve_shape {
    enum struct type_t { flat, ribbon, cylinder };
    enum struct basis_t { bezier, bspline };
    vector<vec3f> P          = {};
    basis_t       basis      = basis_t::bezier;
    int           degree     = 3;
    type_t        type       = type_t::flat;
    vector<vec3f> N          = {};
    float         width0     = 1;
    float         width1     = 1;
    int           splitdepth = 3;
};
struct pbrt_loopsubdiv_shape {
    int           levels  = 3;
    vector<vec3i> indices = {};
    vector<vec3f> P       = {};
};
struct pbrt_nurbs_shape {
    int           nu     = -1;
    int           nv     = -1;
    vector<float> uknots = {};
    vector<float> vknots = {};
    float         u0     = -1;
    float         v0     = -1;
    float         u1     = -1;
    float         v1     = -1;
    vector<vec3f> P      = {};
    vector<float> Pw     = {};
};
struct pbrt_sphere_shape {
    float radius = 1;
    float zmin   = -radius;
    float zmax   = radius;
    float phimax = 360;
};
struct pbrt_disk_shape {
    float height      = 0;
    float radius      = 1;
    float innerradius = 0;
    float phimax      = 360;
};
struct pbrt_cone_shape {
    float radius = 1;
    float height = 1;
    float phimax = 360;
};
struct pbrt_cylinder_shape {
    float radius = 1;
    float zmin   = -1;
    float zmax   = 1;
    float phimax = 360;
};
struct pbrt_hyperboloid_shape {
    vec3f p1     = {0, 0, 0};
    vec3f p2     = {1, 1, 1};
    float phimax = 360;
};
struct pbrt_paraboloid_shape {
    float radius = 1;
    float zmin   = 0;
    float zmax   = 1;
    float phimax = 360;
};
struct pbrt_heightfield_shape {
    int           nu = 0;
    int           nv = 0;
    vector<float> Pz = {};
};
using pbrt_shape = variant<pbrt_trianglemesh_shape, pbrt_plymesh_shape,
    pbrt_curve_shape, pbrt_loopsubdiv_shape, pbrt_nurbs_shape,
    pbrt_sphere_shape, pbrt_disk_shape, pbrt_cone_shape, pbrt_cylinder_shape,
    pbrt_hyperboloid_shape, pbrt_paraboloid_shape, pbrt_heightfield_shape>;

// pbrt lights
struct pbrt_distant_light {
    pbrt_spectrum3f scale = {1, 1, 1};
    pbrt_spectrum3f L     = {1, 1, 1};
    vec3f           from{0, 0, 0};
    vec3f           to = {0, 0, 1};
};
struct pbrt_goniometric_light {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f I       = {1, 1, 1};
    string          mapname = "";
};
struct pbrt_infinite_light {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f L       = {1, 1, 1};
    int             samples = 1;
    string          mapname = "";
};
struct pbrt_point_light {
    pbrt_spectrum3f scale = {1, 1, 1};
    pbrt_spectrum3f I     = {1, 1, 1};
    vec3f           from{0, 0, 0};
};
struct pbrt_projection_light {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f I       = {1, 1, 1};
    float           fov     = 45;
    string          mapname = "";
};
struct pbrt_spot_light {
    pbrt_spectrum3f scale          = {1, 1, 1};
    pbrt_spectrum3f I              = {1, 1, 1};
    vec3f           from           = {0, 0, 0};
    vec3f           to             = {0, 0, 1};
    float           coneangle      = 30;
    float           conedeltaangle = 5;
};
using pbrt_light =
    variant<pbrt_distant_light, pbrt_goniometric_light, pbrt_infinite_light,
        pbrt_point_light, pbrt_projection_light, pbrt_spot_light>;

// pbrt area lights
struct pbrt_none_arealight {};
struct pbrt_diffuse_arealight {
    pbrt_spectrum3f scale    = {1, 1, 1};
    pbrt_spectrum3f L        = {1, 1, 1};
    bool            twosided = false;
    int             samples  = 1;
};
enum struct pbrt_arealight_type { none, diffuse };
struct pbrt_arealight {
    pbrt_arealight_type type = pbrt_arealight_type::none;
    pbrt_none_arealight none = {};
    pbrt_diffuse_arealight diffuse = {};
};

// pbrt mediums
struct pbrt_homogeneous_medium {
    pbrt_spectrum3f sigma_a = {0.0011f, 0.0024f, 0.014f};
    pbrt_spectrum3f sigma_s = {2.55f, 3.21f, 3.77f};
    string          preset  = "";
    float           g       = 0;
    float           scale   = 1;
};
struct pbrt_heterogeneous_medium {
    pbrt_spectrum3f sigma_a = {0.0011f, 0.0024f, 0.014f};
    pbrt_spectrum3f sigma_s = {2.55f, 3.21f, 3.77f};
    string          preset  = "";
    float           g       = 0;
    float           scale   = 1;
    vec3f           p0      = {0, 0, 0};
    vec3f           p1      = {1, 1, 1};
    int             nx      = 1;
    int             ny      = 1;
    int             nz      = 1;
    vector<float>   density = {};
};
enum struct pbrt_medium_type { homogeneous, heterogeneous };
struct pbrt_medium {
    pbrt_medium_type type = pbrt_medium_type::homogeneous;
    pbrt_homogeneous_medium homogeneous = {};
    pbrt_heterogeneous_medium heterogeneous = {};
};

// pbrt medium interface
struct pbrt_mediuminterface {
    string interior = "";
    string exterior = "";
};

// pbrt insstance
struct pbrt_object {
    string name = "";
};

// pbrt stack ctm
struct pbrt_context {
    affine3f transform_start        = identity_affine3f;
    affine3f transform_end          = identity_affine3f;
    string   material               = "";
    string   arealight              = "";
    string   medium_interior        = "";
    string   medium_exterior        = "";
    bool     reverse                = false;
    bool     active_transform_start = true;
    bool     active_transform_end   = true;
    float    last_lookat_distance   = 0;
};

// pbrt callbacks
struct pbrt_callbacks {
    virtual void sampler(const pbrt_sampler& value, const pbrt_context& ctx) {}
    virtual void integrator(const pbrt_integrator& value, const pbrt_context& ctx) {}
    virtual void accelerator(const pbrt_accelerator& value, const pbrt_context& ctx) {}
    virtual void film(const pbrt_film& value, const pbrt_context& ctx) {}
    virtual void filter(const pbrt_filter& value, const pbrt_context& ctx) {}
    virtual void camera(const pbrt_camera& value, const pbrt_context& ctx) {}
    virtual void texture(const pbrt_texture& value, const string& name,
        const pbrt_context& ctx) {}
    virtual void material(const pbrt_material& value, const string& name,
        const pbrt_context& ctx) {}
    virtual void medium(const pbrt_medium& value, const string& name,
        const pbrt_context& ctx) {}
    virtual void shape(const pbrt_shape& value, const pbrt_context& ctx) {}
    virtual void light(const pbrt_light& value, const pbrt_context& ctx) {}
    virtual void arealight(const pbrt_arealight& value, const string& name,
        const pbrt_context& ctx) {}
    virtual void object_instance(const pbrt_object& value, const pbrt_context& ctx) {}
    virtual void begin_object(const pbrt_object& value, const pbrt_context& ctx) {}
    virtual void end_object(const pbrt_object& value, const pbrt_context& ctx) {}
};

// Load pbrt params
struct pbrt_params {
    bool geometry_only = false;
    bool flip_texcoord = true;
};

// Load pbrt scene
void load_pbrt(
    const string& filename, pbrt_callbacks& cb, const pbrt_params& params = {});

}  // namespace yocto

#endif
