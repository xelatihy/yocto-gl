//
// # Yocto/Pbrt: Tiny library for Pbrt parsing
//
// Yocto/Pbrt is a simple pbrt parser that works with callbacks.
// We make no attempt to provide a simple interface for pbrt but just the
// low level parsing code.
//
// Error reporting is done through exceptions using the `pbrtio_error`
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
using std::string_view;
using std::unordered_map;
using std::variant;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// pbrt spectrum as rgb color
template <typename T, int N>
struct spectrum;

// pbrt spectrum as rgb color
template <typename T>
struct spectrum<T, 3> {
    union {
        struct {
            T x, y, z;
        };
        T elems[3];
    };

    constexpr spectrum() : x{0}, y{0}, z{0} {}
    constexpr spectrum(T x, T y, T z) : x{x}, y{y}, z{z} {}
    constexpr explicit spectrum(T v) : x{v}, y{v}, z{v} {}
    constexpr explicit operator vec<T, 3>() const { return {x, y, z}; };

    constexpr T&       operator[](int i) { return elems[i]; }
    constexpr const T& operator[](int i) const { return elems[i]; }
};

// typedefs
using spectrum3f = spectrum<float, 3>;

// pbrt cameras
struct pbrt_camera_perspective {
    float  fov              = 90;
    float  frameaspectratio = -1;  // or computed from film
    float  lensradius       = 0;
    float  focaldistance    = 1e30;
    bbox2f screenwindow     = {{-1, -1}, {1, 1}};
    float  shutteropen      = 0;
    float  shutterclose     = 1;
};
struct pbrt_camera_orthographic {
    float  frameaspectratio = -1;  // or computed from film
    float  lensradius       = 0;
    float  focaldistance    = 1e30;
    bbox2f screenwindow     = {{-1, -1}, {1, 1}};
    float  shutteropen      = 0;
    float  shutterclose     = 1;
};
struct pbrt_camera_environment {
    float shutteropen  = 0;
    float shutterclose = 1;
};
struct pbrt_camera_realistic {
    string lensfile         = "";
    float  aperturediameter = 1;
    float  focusdistance    = 10;
    bool   simpleweighting  = true;
    float  shutteropen      = 0;
    float  shutterclose     = 1;
};
using pbrt_camera = variant<pbrt_camera_perspective, pbrt_camera_orthographic,
    pbrt_camera_environment, pbrt_camera_realistic>;

// pbrt samplers
struct pbrt_sampler_random {
    int pixelsamples = 16;
};
struct pbrt_sampler_halton {
    int pixelsamples = 16;
};
struct pbrt_sampler_sobol {
    int pixelsamples = 16;
};
struct pbrt_sampler_zerotwosequence {
    int pixelsamples = 16;
};
struct pbrt_sampler_maxmindist {
    int pixelsamples = 16;
};
struct pbrt_sampler_stratified {
    bool jitter   = true;
    int  xsamples = 2;
    int  ysamples = 2;
};
using pbrt_sampler = variant<pbrt_sampler_random, pbrt_sampler_halton,
    pbrt_sampler_sobol, pbrt_sampler_zerotwosequence, pbrt_sampler_maxmindist,
    pbrt_sampler_stratified>;

// pbrt film
struct pbrt_film_image {
    int    xresolution        = 640;
    int    yresolution        = 480;
    bbox2f cropwindow         = {{0, 0}, {1, 1}};
    float  scale              = 1;
    float  maxsampleluminance = type_max<float>;
    float  diagonal           = 35;
    string filename           = "pbrt.exr";
};
using pbrt_film = variant<pbrt_film_image>;

// pbrt filters
struct pbrt_filter_box {
    float xwidth = 0.5f;
    float ywidth = 0.5f;
};
struct pbrt_filter_gaussian {
    float xwidth = 2;
    float ywidth = 2;
    float alpha  = 2;
};
struct pbrt_filter_mitchell {
    float xwidth = 2;
    float ywidth = 2;
    float B      = 1.0f / 3.0f;
    float C      = 1.0f / 3.0f;
};
struct pbrt_filter_sinc {
    float xwidth = 4;
    float ywidth = 4;
    float tau    = 3;
};
struct pbrt_filter_triangle {
    float xwidth = 2;
    float ywidth = 2;
};
using pbrt_filter = variant<pbrt_filter_box, pbrt_filter_gaussian,
    pbrt_filter_mitchell, pbrt_filter_sinc, pbrt_filter_triangle>;

// pbrt integrators
struct pbrt_integrator_path {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    float  rrthreshold = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
};
struct pbrt_integrator_volpath {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    float  rrthreshold = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
};
struct pbrt_integrator_bdpt {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::power;
    bool                  visualizestrategies = false;
    bool                  visualizeweights    = false;
};
struct pbrt_integrator_directlighting {
    enum struct strategy_t { all, one };
    strategy_t strategy    = strategy_t::all;
    int        maxdepth    = 5;
    bbox2i     pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
};
struct pbrt_integrator_mlt {
    int    maxdepth             = 5;
    bbox2i pixelbounds          = {{0, 0}, {type_max<int>, type_max<int>}};
    int    bootstrapsamples     = 100000;
    int    chains               = 1000;
    int    mutationsperpixel    = 100;
    float  largestepprobability = 0.3;
    float  sigma                = 0.01;
};
struct pbrt_integrator_sppm {
    int    maxdepth            = 5;
    bbox2i pixelbounds         = {{0, 0}, {type_max<int>, type_max<int>}};
    int    iterations          = 64;
    int    photonsperiteration = -1;
    int    imagewritefrequency = pow2(31);
    float  radius              = 5;
};
struct pbrt_integrator_whitted {
    int    maxdepth    = 5;
    bbox2i pixelbounds = {{0, 0}, {type_max<int>, type_max<int>}};
};
using pbrt_integrator = variant<pbrt_integrator_path, pbrt_integrator_volpath,
    pbrt_integrator_bdpt, pbrt_integrator_directlighting, pbrt_integrator_mlt,
    pbrt_integrator_sppm, pbrt_integrator_whitted>;

// pbrt accellerators
struct pbrt_accelerator_bvh {
    enum struct splitmethod_t { sah, equal, middle, hlbvh };
    int           maxnodeprims = 4;
    splitmethod_t splitmethod  = splitmethod_t::sah;
};
struct pbrt_accelerator_kdtree {
    int   intersectcost = 80;
    int   traversalcost = 1;
    float emptybonus    = 0.2;
    int   maxprims      = 1;
    int   maxdepth      = -1;
};
using pbrt_accelerator = variant<pbrt_accelerator_bvh, pbrt_accelerator_kdtree>;

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
struct pbrt_textured<spectrum3f> {
    spectrum3f value   = {0, 0, 0};
    string     texture = "";
    pbrt_textured() : value{0, 0, 0}, texture{} {}
    pbrt_textured(float x, float y, float z) : value{x, y, z}, texture{} {}
};

// pbrt textures
struct pbrt_texture_constant {
    pbrt_textured<spectrum3f> value = {1, 1, 1};
};
struct pbrt_texture_bilerp {
    pbrt_textured<spectrum3f> v00 = {0, 0, 0};
    pbrt_textured<spectrum3f> v01 = {1, 1, 1};
    pbrt_textured<spectrum3f> v10 = {0, 0, 0};
    pbrt_textured<spectrum3f> v11 = {1, 1, 1};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_checkerboard {
    enum struct aamode_type { closedform, none };
    int                       dimension = 2;
    pbrt_textured<spectrum3f> tex1      = {1, 1, 1};
    pbrt_textured<spectrum3f> tex2      = {0, 0, 0};
    aamode_type               aamode    = aamode_type::closedform;
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_dots {
    pbrt_textured<spectrum3f> inside  = {1, 1, 1};
    pbrt_textured<spectrum3f> outside = {0, 0, 0};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_fbm {
    int   octaves   = 8;
    float roughness = 0.5;
};
struct pbrt_texture_imagemap {
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
struct pbrt_texture_marble {
    int   octaves   = 8;
    float roughness = 0.5f;
    float scale     = 1;
    float variation = 0.2f;
};
struct pbrt_texture_mix {
    pbrt_textured<spectrum3f> tex1   = {1, 1, 1};
    pbrt_textured<spectrum3f> tex2   = {1, 1, 1};
    pbrt_textured<float>      amount = 0.5f;
};
struct pbrt_texture_scale {
    pbrt_textured<spectrum3f> tex1 = {1, 1, 1};
    pbrt_textured<spectrum3f> tex2 = {1, 1, 1};
};
struct pbrt_texture_uv {
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_windy {
    // TODO: missing parameters
};
struct pbrt_texture_wrinkled {
    int   octaves   = 8;
    float roughness = 0.5;
};
using pbrt_texture = variant<pbrt_texture_constant, pbrt_texture_bilerp,
    pbrt_texture_checkerboard, pbrt_texture_dots, pbrt_texture_fbm,
    pbrt_texture_imagemap, pbrt_texture_marble, pbrt_texture_mix,
    pbrt_texture_scale, pbrt_texture_uv, pbrt_texture_windy,
    pbrt_texture_wrinkled>;

// pbrt materials
struct pbrt_material_matte {
    pbrt_textured<spectrum3f> Kd      = {0.5f, 0.5f, 0.5f};
    pbrt_textured<float>      sigma   = 0;
    pbrt_textured<float>      bumpmap = 0;
};
struct pbrt_material_mirror {
    pbrt_textured<spectrum3f> Kr      = {0.9f, 0.9f, 0.9f};
    pbrt_textured<float>      bumpmap = 0;
};
struct pbrt_material_plastic {
    pbrt_textured<spectrum3f> Kd             = {0.25f, 0.25f, 0.25f};
    pbrt_textured<spectrum3f> Ks             = {0.25f, 0.25f, 0.25f};
    pbrt_textured<float>      uroughness     = 0.1f;
    pbrt_textured<float>      vroughness     = 0.1f;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_metal {
    pbrt_textured<spectrum3f> eta = {
        0.2004376970f, 0.9240334304f, 1.1022119527f};
    pbrt_textured<spectrum3f> k = {3.9129485033f, 2.4528477015f, 2.1421879552f};
    pbrt_textured<float>      uroughness     = 0.01;
    pbrt_textured<float>      vroughness     = 0.01;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_glass {
    pbrt_textured<spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>      eta            = 1;
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_translucent {
    pbrt_textured<spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<spectrum3f> reflect        = {0, 0, 0};
    pbrt_textured<spectrum3f> transmit       = {0, 0, 0};
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_uber {
    pbrt_textured<spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<spectrum3f> Kr             = {0, 0, 0};
    pbrt_textured<spectrum3f> Kt             = {0, 0, 0};
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    pbrt_textured<float>      eta            = 1;
    pbrt_textured<spectrum3f> opacity        = {1, 1, 1};
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_disney {
    pbrt_textured<spectrum3f> color           = {0.5f, 0.5f, 0.5f};
    pbrt_textured<float>      anisotropic     = 0;
    pbrt_textured<float>      clearcoat       = 0;
    pbrt_textured<float>      clearcoatgloss  = 1;
    pbrt_textured<float>      eta             = 1.5f;
    pbrt_textured<float>      metallic        = 0;
    pbrt_textured<float>      uroughness      = 0;
    pbrt_textured<float>      vroughness      = 0;
    pbrt_textured<spectrum3f> scatterdistance = {0, 0, 0};
    pbrt_textured<float>      sheen           = 0;
    pbrt_textured<float>      sheentint       = 0.5;
    pbrt_textured<float>      spectrans       = 0;
    pbrt_textured<float>      speculartint    = 0;
    bool                      thin            = false;
    pbrt_textured<spectrum3f> difftrans       = {1, 1, 1};
    pbrt_textured<spectrum3f> flatness        = {0, 0, 0};
    bool                      remaproughness  = true;
    pbrt_textured<float>      bumpmap         = 0;
};
struct pbrt_material_fourier {
    string               bsdffile = "";
    pbrt_textured<float> bumpmap  = 0;
};
struct pbrt_material_hair {
    pbrt_textured<spectrum3f> color       = {0, 0, 0};  // TODO: missing default
    pbrt_textured<spectrum3f> sigma_a     = {0, 0, 0};  // TODO: missing default
    pbrt_textured<float>      eumelanin   = 0;          // TODO: missing default
    pbrt_textured<float>      pheomelanin = 0;          // TODO: missing default
    pbrt_textured<float>      eta         = 1.55f;
    pbrt_textured<float>      beta_m      = 0.3f;
    pbrt_textured<float>      beta_n      = 0.3f;
    pbrt_textured<float>      alpha       = 2;
    pbrt_textured<float>      bumpmap     = 0;
};
struct pbrt_material_kdsubsurface {
    pbrt_textured<spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<spectrum3f> mfp            = {1, 1, 1};
    pbrt_textured<float>      eta            = 1;
    pbrt_textured<spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_mix {
    pbrt_textured<spectrum3f> amount         = {0, 0, 0};
    string                    namedmaterial1 = "";
    string                    namedmaterial2 = "";
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_substrate {
    pbrt_textured<spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
struct pbrt_material_subsurface {
    string                    name           = "";
    pbrt_textured<spectrum3f> sigma_a        = {.0011, .0024, .014};
    pbrt_textured<spectrum3f> sigma_prime_s  = {2.55, 3.12, 3.77};
    float                     scale          = 1;
    pbrt_textured<float>      eta            = 1;
    pbrt_textured<spectrum3f> Kr             = {1, 1, 1};
    pbrt_textured<spectrum3f> Kt             = {1, 1, 1};
    pbrt_textured<float>      uroughness     = 0;
    pbrt_textured<float>      vroughness     = 0;
    bool                      remaproughness = true;
    pbrt_textured<float>      bumpmap        = 0;
};
using pbrt_material = variant<pbrt_material_matte, pbrt_material_mirror,
    pbrt_material_plastic, pbrt_material_metal, pbrt_material_glass,
    pbrt_material_translucent, pbrt_material_uber, pbrt_material_disney,
    pbrt_material_fourier, pbrt_material_hair, pbrt_material_kdsubsurface,
    pbrt_material_mix, pbrt_material_substrate, pbrt_material_subsurface>;

// pbrt shapes
struct pbrt_shape_trianglemesh {
    vector<vec3i>        indices     = {};
    vector<vec3f>        P           = {};
    vector<vec3f>        N           = {};
    vector<vec3f>        S           = {};
    vector<vec2f>        uv          = {};
    pbrt_textured<float> alpha       = 1;
    pbrt_textured<float> shadowalpha = 1;
};
struct pbrt_shape_plymesh {
    string               filename    = {};
    pbrt_textured<float> alpha       = 1;
    pbrt_textured<float> shadowalpha = 1;
};
struct pbrt_shape_curve {
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
struct pbrt_shape_loopsubdiv {
    int           levels  = 3;
    vector<vec3i> indices = {};
    vector<vec3f> P       = {};
};
struct pbrt_shape_nurbs {
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
struct pbrt_shape_sphere {
    float radius = 1;
    float zmin   = -radius;
    float zmax   = radius;
    float phimax = 360;
};
struct pbrt_shape_disk {
    float height      = 0;
    float radius      = 1;
    float innerradius = 0;
    float phimax      = 360;
};
struct pbrt_shape_cone {
    float radius = 1;
    float height = 1;
    float phimax = 360;
};
struct pbrt_shape_cylinder {
    float radius = 1;
    float zmin   = -1;
    float zmax   = 1;
    float phimax = 360;
};
struct pbrt_shape_hyperboloid {
    vec3f p1     = {0, 0, 0};
    vec3f p2     = {1, 1, 1};
    float phimax = 360;
};
struct pbrt_shape_paraboloid {
    float radius = 1;
    float zmin   = 0;
    float zmax   = 1;
    float phimax = 360;
};
struct pbrt_shape_heightfield {
    int           nu = 0;
    int           nv = 0;
    vector<float> Pz = {};
};
using pbrt_shape = variant<pbrt_shape_trianglemesh, pbrt_shape_plymesh,
    pbrt_shape_curve, pbrt_shape_loopsubdiv, pbrt_shape_nurbs,
    pbrt_shape_sphere, pbrt_shape_disk, pbrt_shape_cone, pbrt_shape_cylinder,
    pbrt_shape_hyperboloid, pbrt_shape_paraboloid, pbrt_shape_heightfield>;

// pbrt lights
struct pbrt_light_distant {
    spectrum3f scale = {1, 1, 1};
    spectrum3f L     = {1, 1, 1};
    vec3f      from{0, 0, 0};
    vec3f      to = {0, 0, 1};
};
struct pbrt_light_goniometric {
    spectrum3f scale   = {1, 1, 1};
    spectrum3f I       = {1, 1, 1};
    string     mapname = "";
};
struct pbrt_light_infinite {
    spectrum3f scale   = {1, 1, 1};
    spectrum3f L       = {1, 1, 1};
    int        samples = 1;
    string     mapname = "";
};
struct pbrt_light_point {
    spectrum3f scale = {1, 1, 1};
    spectrum3f I     = {1, 1, 1};
    vec3f      from{0, 0, 0};
};
struct pbrt_light_projection {
    spectrum3f scale   = {1, 1, 1};
    spectrum3f I       = {1, 1, 1};
    float      fov     = 45;
    string     mapname = "";
};
struct pbrt_light_spot {
    spectrum3f scale          = {1, 1, 1};
    spectrum3f I              = {1, 1, 1};
    vec3f      from           = {0, 0, 0};
    vec3f      to             = {0, 0, 1};
    float      coneangle      = 30;
    float      conedeltaangle = 5;
};
using pbrt_light =
    variant<pbrt_light_distant, pbrt_light_goniometric, pbrt_light_infinite,
        pbrt_light_point, pbrt_light_projection, pbrt_light_spot>;

// pbrt area lights
struct pbrt_arealight_none {};
struct pbrt_arealight_diffuse {
    spectrum3f scale    = {1, 1, 1};
    spectrum3f L        = {1, 1, 1};
    bool       twosided = false;
    int        samples  = 1;
};
using pbrt_arealight = variant<pbrt_arealight_none, pbrt_arealight_diffuse>;

// pbrt mediums
struct pbrt_medium_homogeneous {
    spectrum3f sigma_a = {0.0011f, 0.0024f, 0.014f};
    spectrum3f sigma_s = {2.55f, 3.21f, 3.77f};
    string     preset  = "";
    float      g       = 0;
    float      scale   = 1;
};
struct pbrt_medium_heterogeneous {
    spectrum3f    sigma_a = {0.0011f, 0.0024f, 0.014f};
    spectrum3f    sigma_s = {2.55f, 3.21f, 3.77f};
    string        preset  = "";
    float         g       = 0;
    float         scale   = 1;
    vec3f         p0      = {0, 0, 0};
    vec3f         p1      = {1, 1, 1};
    int           nx      = 1;
    int           ny      = 1;
    int           nz      = 1;
    vector<float> density = {};
};
using pbrt_medium = variant<pbrt_medium_homogeneous, pbrt_medium_heterogeneous>;

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
    affine3f frame           = identity_affine3f;
    string   material        = "";
    string   arealight       = "";
    string   medium_interior = "";
    string   medium_exterior = "";
    bool     reverse         = false;
};

// pbrt callbacks
struct pbrt_callbacks {
    void sampler(const pbrt_sampler& value, const pbrt_context& ctx) {}
    void integrator(const pbrt_integrator& value, const pbrt_context& ctx) {}
    void film(const pbrt_film& value, const pbrt_context& ctx) {}
    void filter(const pbrt_filter& value, const pbrt_context& ctx) {}
    void camera(const pbrt_camera& value, const pbrt_context& ctx) {}
    void texture(const pbrt_texture& value, const string& name,
        const pbrt_context& ctx) {}
    void material(const pbrt_material& value, const string& name,
        const pbrt_context& ctx) {}
    void medium(const pbrt_medium& value, const string& name,
        const pbrt_context& ctx) {}
    void shape(const pbrt_shape& value, const pbrt_context& ctx) {}
    void light(const pbrt_light& value, const pbrt_context& ctx) {}
    void arealight(const pbrt_arealight& value, const string& name,
        const pbrt_context& ctx) {}
    void object_instance(const pbrt_object& value, const pbrt_context& ctx) {}
    void begin_object(const pbrt_object& value, const pbrt_context& ctx) {}
    void end_object(const pbrt_object& value, const pbrt_context& ctx) {}
};

// Load obj options
struct load_pbrt_options {
    bool geometry_only = false;
    bool flip_texcoord = true;
};

// Load obj scene
template <typename Callbacks>
inline void load_pbrt(const string& filename, Callbacks& cb,
    const load_pbrt_options& options = {});

// objio error
struct pbrtio_error : runtime_error {
    explicit pbrtio_error(const char* msg) : runtime_error{msg} {}
    explicit pbrtio_error(const std::string& msg) : runtime_error{msg} {}
};

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF LOW LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Token stream
struct pbrt_token_stream {
    string      buffer;
    string_view str;
    string_view saved;
};

// skip white space or comment
static inline void skip_whitespace_or_comment(pbrt_token_stream& stream) {
    auto& str = stream.str;
    if (str.empty()) return;
    while (!str.empty() && (std::isspace(str.front()) || str.front() == '#' ||
                               str.front() == ',')) {
        if (str.front() == '#') {
            auto pos = str.find('\n');
            if (pos != string_view::npos) {
                str.remove_prefix(pos);
            } else {
                str.remove_prefix(str.length());
            }
        } else {
            auto pos = str.find_first_not_of(" \t\n\r,");
            if (pos == string_view::npos) {
                str.remove_prefix(str.length());
            } else {
                str.remove_prefix(pos);
            }
        }
    }
}

// parse a quoted string
static inline void parse_value(pbrt_token_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.front() != '"') {
        throw pbrtio_error("bad string");
    }
    str.remove_prefix(1);
    auto pos = str.find('"');
    if (pos == string_view::npos) {
        throw pbrtio_error("bad string");
    }
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
}

// parse a quoted string
static inline void parse_command(pbrt_token_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (!std::isalpha((int)str.front())) {
        throw pbrtio_error("bad command");
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
static inline void parse_value(pbrt_token_stream& stream, float& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw pbrtio_error("number expected");
    auto next = (char*)nullptr;
    value     = strtof(str.data(), &next);
    if (str.data() == next) throw pbrtio_error("number expected");
    str.remove_prefix(next - str.data());
}

// parse a number
static inline void parse_value(pbrt_token_stream& stream, int& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw pbrtio_error("number expected");
    auto next = (char*)nullptr;
    value     = strtol(str.data(), &next, 10);
    if (str.data() == next) throw pbrtio_error("number expected");
    str.remove_prefix(next - str.data());
}
static inline void parse_value(pbrt_token_stream& stream, bool& value) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    if (value_name == "true") {
        value = true;
    } else if (value_name == "false") {
        value = false;
    } else {
        throw pbrtio_error("expected boolean");
    }
}
template <typename T>
static inline void parse_value(pbrt_token_stream& stream, T& value,
    unordered_map<string, T>& value_names) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw pbrtio_error("expected enum value");
    }
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_bilerp::mapping_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_bilerp::mapping_type>{
            {"uv", pbrt_texture_bilerp::mapping_type::uv},
            {"spherical", pbrt_texture_bilerp::mapping_type::spherical},
            {"cylindrical", pbrt_texture_bilerp::mapping_type::cylindrical},
            {"planar", pbrt_texture_bilerp::mapping_type::planar},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_checkerboard::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_dots::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_imagemap::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_uv::mapping_type& value) {
    return parse_value(stream, (pbrt_texture_bilerp::mapping_type&)value);
}

static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_checkerboard::aamode_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_checkerboard::aamode_type>{
            {"closedform", pbrt_texture_checkerboard::aamode_type::closedform},
            {"none", pbrt_texture_checkerboard::aamode_type::none},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_texture_imagemap::wrap_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_texture_imagemap::wrap_type>{
            {"repeat", pbrt_texture_imagemap::wrap_type::repeat},
            {"clamp", pbrt_texture_imagemap::wrap_type::clamp},
            {"black", pbrt_texture_imagemap::wrap_type::black},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_shape_curve::basis_t& value) {
    static auto value_names = unordered_map<string, pbrt_shape_curve::basis_t>{
        {"bezier", pbrt_shape_curve::basis_t::bezier},
        {"bspline", pbrt_shape_curve::basis_t::bspline},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_token_stream& stream, pbrt_shape_curve::type_t& value) {
    static auto value_names = unordered_map<string, pbrt_shape_curve::type_t>{
        {"flat", pbrt_shape_curve::type_t::flat},
        {"cylinder", pbrt_shape_curve::type_t::cylinder},
        {"ribbon", pbrt_shape_curve::type_t::ribbon},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_path::lightsamplestrategy_t&  value) {
    static auto value_names =
        unordered_map<string, pbrt_integrator_path::lightsamplestrategy_t>{
            {"power", pbrt_integrator_path::lightsamplestrategy_t::power},
            {"spatial", pbrt_integrator_path::lightsamplestrategy_t::spatial},
            {"uniform", pbrt_integrator_path::lightsamplestrategy_t::uniform},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(pbrt_token_stream&   stream,
    pbrt_integrator_volpath::lightsamplestrategy_t& value) {
    return parse_value(
        stream, (pbrt_integrator_path::lightsamplestrategy_t&)value);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_bdpt::lightsamplestrategy_t&  value) {
    return parse_value(
        stream, (pbrt_integrator_path::lightsamplestrategy_t&)value);
}
static inline void parse_value(pbrt_token_stream& stream,
    pbrt_integrator_directlighting::strategy_t&   value) {
    static auto value_names =
        unordered_map<string, pbrt_integrator_directlighting::strategy_t>{
            {"all", pbrt_integrator_directlighting::strategy_t::all},
            {"one", pbrt_integrator_directlighting::strategy_t::one},
        };
    return parse_value(stream, value, value_names);
}

// parse a vec type
template <typename T, int N>
static inline void parse_value(pbrt_token_stream& stream, vec<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}
template <typename T, int N, int M>
static inline void parse_value(pbrt_token_stream& stream, mat<T, N, M>& value) {
    for (auto i = 0; i < M; i++) parse_value(stream, value[i]);
}
template <typename T>
static inline void parse_value(pbrt_token_stream& stream, bbox<T, 2>& value) {
    parse_value(stream, value[0][0]);
    parse_value(stream, value[1][0]);
    parse_value(stream, value[0][1]);
    parse_value(stream, value[1][1]);
}
template <typename T, int N>
static inline void parse_value(
    pbrt_token_stream& stream, spectrum<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}

// Check next
static inline bool is_empty(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return stream.str.empty();
}
static inline bool is_string(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '"';
}
static inline bool is_open_bracket(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '[';
}
static inline bool is_close_bracket(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == ']';
}
static inline bool is_param(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    return is_string(stream);
}

// parse a quoted string
static inline void parse_nametype(
    pbrt_token_stream& stream, string& name, string& type) {
    auto value = ""s;
    parse_value(stream, value);
    auto str  = string_view{value};
    auto pos1 = str.find(' ');
    if (pos1 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    type = string(str.substr(0, pos1));
    str.remove_prefix(pos1);
    auto pos2 = str.find_first_not_of(' ');
    if (pos2 == string_view::npos) {
        throw pbrtio_error("bad type " + value);
    }
    str.remove_prefix(pos2);
    name = string(str);
}

static inline void skip_open_bracket(pbrt_token_stream& stream) {
    if (!is_open_bracket(stream)) throw pbrtio_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}
static inline void skip_close_bracket(pbrt_token_stream& stream) {
    if (!is_close_bracket(stream)) throw pbrtio_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}

template <typename T>
static inline void parse_param(pbrt_token_stream& stream, T& value) {
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}

template <typename T>
static inline void parse_param(pbrt_token_stream& stream, vector<T>& values) {
    skip_open_bracket(stream);
    values.clear();
    while (!is_close_bracket(stream)) {
        values.push_back({});
        parse_value(stream, values.back());
    }
    skip_close_bracket(stream);
}

template <typename T>
static inline bool is_type_compatible(const string& type) {
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
    } else if constexpr (std::is_same<T, spectrum3f>::value) {
        return type == "rgb" || type == "spectrum" || type == "blackbody";
    } else if constexpr (std::is_same<T, vec3i>::value) {
        return type == "integer";
    } else if constexpr (std::is_same<T, bbox2i>::value) {
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
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, T& value) {
    if (!is_type_compatible<T>(type)) {
        throw pbrtio_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, spectrum<T, 3>& value) {
    bool verbose = false;
    if (type == "rgb") {
        parse_param(stream, value);
    } else if (type == "color") {
        parse_param(stream, value);
    } else if (type == "float") {
        auto valuef = 0.0f;
        parse_param(stream, valuef);
        value = {valuef, valuef, valuef};
    } else if (type == "blackbody") {
        auto blackbody = zero2f;
        parse_param(stream, blackbody);
        (vec3f&)value = blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "spectrum" && is_string(stream)) {
        if (verbose) printf("spectrum  not well supported\n");
        auto filename = ""s;
        parse_param(stream, filename);
        value = {1, 0, 0};
    } else if (type == "spectrum" && !is_string(stream)) {
        if (verbose) printf("spectrum  not well supported\n");
        auto values = vector<float>{};
        parse_param(stream, values);
        value = {1, 0, 0};
    } else {
        throw pbrtio_error("unsupported spectrum type");
    }
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, vector<T>& value) {
    if (!is_type_compatible<T>(type)) {
        throw pbrtio_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

template <typename T>
static inline void parse_param(
    pbrt_token_stream& stream, const string& type, pbrt_textured<T>& value) {
    if (type == "texture") {
        parse_param(stream, value.texture);
    } else {
        parse_param(stream, type, value.value);
    }
}

static inline void skip_value(pbrt_token_stream& stream) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.front() == '"') {
        str.remove_prefix(1);
        str.remove_prefix(str.find('"') + 1);
    } else {
        str.remove_prefix(str.find_first_of(" \n\t\r],\""));
    }
    skip_whitespace_or_comment(stream);
}

static inline void skip_param(pbrt_token_stream& stream) {
    if (is_open_bracket(stream)) {
        skip_open_bracket(stream);
        while (!is_close_bracket(stream)) skip_value(stream);
        skip_close_bracket(stream);
    } else {
        skip_value(stream);
    }
}

static inline void save_stream_position(pbrt_token_stream& stream) {
    stream.saved = stream.str;
}
static inline void restore_stream_position(pbrt_token_stream& stream) {
    stream.str = stream.saved;
}

// Load a token stream
static inline void load_token_stream(
    const string& filename, pbrt_token_stream& stream) {
    load_text(filename, stream.buffer);
    stream.str = stream.buffer;
}

// operations on token stacks
static inline void init_token_stream(vector<pbrt_token_stream>& stream) {
    stream.reserve(100);
}
static inline void load_token_stream(
    const string& filename, vector<pbrt_token_stream>& stream) {
    stream.emplace_back();
    load_token_stream(filename, stream.back());
}

// Skip whitespace
static inline void skip_whitespace_or_comment_to_next_file(
    vector<pbrt_token_stream>& stream) {
    if (stream.empty()) return;
    while (!stream.empty()) {
        skip_whitespace_or_comment(stream.back());
        if (is_empty(stream.back())) {
            stream.pop_back();
        } else {
            break;
        }
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Parse Integrator
static inline void parse_pbrt_integrator(
    pbrt_token_stream& stream, const string& type, pbrt_integrator& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "path") {
        auto tvalue = pbrt_integrator_path{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "rrthreshold") {
                parse_param(stream, ptype, tvalue.rrthreshold);
            } else if (pname == "lightsamplestrategy") {
                parse_param(stream, ptype, tvalue.lightsamplestrategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
            // parse_optional_param(stream, "lightsamplestrategy",
            // tvalue.lightsamplestrategy); // TODO: enums
        }
        value = tvalue;
    } else if (type == "volpath") {
        auto tvalue = pbrt_integrator_volpath{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "rrthreshold") {
                parse_param(stream, ptype, tvalue.rrthreshold);
            } else if (pname == "lightsamplestrategy") {
                parse_param(stream, ptype, tvalue.lightsamplestrategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "directlighting") {
        auto tvalue = pbrt_integrator_directlighting{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "strategy") {
                parse_param(stream, ptype, tvalue.strategy);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bdpt") {
        auto tvalue = pbrt_integrator_bdpt{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "lightsamplestrategy") {
                parse_param(stream, ptype, tvalue.lightsamplestrategy);
            } else if (pname == "visualizestrategies") {
                parse_param(stream, ptype, tvalue.visualizestrategies);
            } else if (pname == "visualizeweights") {
                parse_param(stream, ptype, tvalue.visualizeweights);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mlt") {
        auto tvalue = pbrt_integrator_mlt{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "bootstrapsamples") {
                parse_param(stream, ptype, tvalue.bootstrapsamples);
            } else if (pname == "chains") {
                parse_param(stream, ptype, tvalue.chains);
            } else if (pname == "mutationsperpixel") {
                parse_param(stream, ptype, tvalue.mutationsperpixel);
            } else if (pname == "largestepprobability") {
                parse_param(stream, ptype, tvalue.largestepprobability);
            } else if (pname == "sigma") {
                parse_param(stream, ptype, tvalue.sigma);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sppm") {
        auto tvalue = pbrt_integrator_sppm{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "iterations") {
                parse_param(stream, ptype, tvalue.iterations);
            } else if (pname == "photonsperiteration") {
                parse_param(stream, ptype, tvalue.photonsperiteration);
            } else if (pname == "imagewritefrequency") {
                parse_param(stream, ptype, tvalue.imagewritefrequency);
            } else if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "whitted") {
        auto tvalue = pbrt_integrator_whitted{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Integrator " + type);
    }
}

// Parse Sampler
static inline void parse_pbrt_sampler(
    pbrt_token_stream& stream, const string& type, pbrt_sampler& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "random") {
        auto tvalue = pbrt_sampler_random{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "halton") {
        auto tvalue = pbrt_sampler_halton{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sobol") {
        auto tvalue = pbrt_sampler_sobol{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "02sequence") {
        auto tvalue = pbrt_sampler_zerotwosequence{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "lowdiscrepancy") {
        auto tvalue = pbrt_sampler_zerotwosequence{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "maxmindist") {
        auto tvalue = pbrt_sampler_maxmindist{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "stratified") {
        auto tvalue = pbrt_sampler_stratified{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xsamples") {
                parse_param(stream, ptype, tvalue.xsamples);
            } else if (pname == "ysamples") {
                parse_param(stream, ptype, tvalue.ysamples);
            } else if (pname == "jitter") {
                parse_param(stream, ptype, tvalue.jitter);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Sampler " + type);
    }
}

// Parse Filter
static inline void parse_pbrt_filter(
    pbrt_token_stream& stream, const string& type, pbrt_filter& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "box") {
        auto tvalue = pbrt_filter_box{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "gaussian") {
        auto tvalue = pbrt_filter_gaussian{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else if (pname == "alpha") {
                parse_param(stream, ptype, tvalue.alpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mitchell") {
        auto tvalue = pbrt_filter_mitchell{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else if (pname == "B") {
                parse_param(stream, ptype, tvalue.B);
            } else if (pname == "C") {
                parse_param(stream, ptype, tvalue.C);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sinc") {
        auto tvalue = pbrt_filter_sinc{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else if (pname == "tau") {
                parse_param(stream, ptype, tvalue.tau);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "triangle") {
        auto tvalue = pbrt_filter_triangle{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown PixelFilter " + type);
    }
}

// Parse Filter
static inline void parse_pbrt_film(
    pbrt_token_stream& stream, const string& type, pbrt_film& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "image") {
        auto tvalue = pbrt_film_image{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xresolution") {
                parse_param(stream, ptype, tvalue.xresolution);
            } else if (pname == "yresolution") {
                parse_param(stream, ptype, tvalue.yresolution);
            } else if (pname == "yresolution") {
                parse_param(stream, ptype, tvalue.yresolution);
            } else if (pname == "cropwindow") {
                parse_param(stream, ptype, tvalue.cropwindow);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "maxsampleluminance") {
                parse_param(stream, ptype, tvalue.maxsampleluminance);
            } else if (pname == "diagonal") {
                parse_param(stream, ptype, tvalue.diagonal);
            } else if (pname == "filename") {
                parse_param(stream, ptype, tvalue.filename);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Camera
static inline void parse_pbrt_camera(
    pbrt_token_stream& stream, const string& type, pbrt_camera& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "perspective") {
        auto tvalue = pbrt_camera_perspective{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "fov") {
                parse_param(stream, ptype, tvalue.fov);
            } else if (pname == "frameaspectratio") {
                parse_param(stream, ptype, tvalue.frameaspectratio);
            } else if (pname == "lensradius") {
                parse_param(stream, ptype, tvalue.lensradius);
            } else if (pname == "focaldistance") {
                parse_param(stream, ptype, tvalue.focaldistance);
            } else if (pname == "screenwindow") {
                parse_param(stream, ptype, tvalue.screenwindow);
            } else if (pname == "shutteropen") {
                parse_param(stream, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(stream, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "orthographic") {
        auto tvalue = pbrt_camera_orthographic{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "frameaspectratio") {
                parse_param(stream, ptype, tvalue.frameaspectratio);
            } else if (pname == "lensradius") {
                parse_param(stream, ptype, tvalue.lensradius);
            } else if (pname == "focaldistance") {
                parse_param(stream, ptype, tvalue.focaldistance);
            } else if (pname == "screenwindow") {
                parse_param(stream, ptype, tvalue.screenwindow);
            } else if (pname == "shutteropen") {
                parse_param(stream, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(stream, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "environment") {
        auto tvalue = pbrt_camera_environment{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "shutteropen") {
                parse_param(stream, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(stream, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "realistic") {
        auto tvalue = pbrt_camera_realistic{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "lensfile") {
                parse_param(stream, ptype, tvalue.lensfile);
            } else if (pname == "aperturediameter") {
                parse_param(stream, ptype, tvalue.aperturediameter);
            } else if (pname == "focusdistance") {
                parse_param(stream, ptype, tvalue.focusdistance);
            } else if (pname == "simpleweighting") {
                parse_param(stream, ptype, tvalue.simpleweighting);
            } else if (pname == "shutteropen") {
                parse_param(stream, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(stream, ptype, tvalue.shutterclose);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Texture
static inline void parse_pbrt_texture(
    pbrt_token_stream& stream, const string& type, pbrt_texture& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "constant") {
        auto tvalue = pbrt_texture_constant{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "value") {
                parse_param(stream, ptype, tvalue.value);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bilerp") {
        auto tvalue = pbrt_texture_bilerp{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "v00") {
                parse_param(stream, ptype, tvalue.v00);
            } else if (pname == "v01") {
                parse_param(stream, ptype, tvalue.v01);
            } else if (pname == "v10") {
                parse_param(stream, ptype, tvalue.v10);
            } else if (pname == "v11") {
                parse_param(stream, ptype, tvalue.v11);
            } else if (pname == "mapping") {
                parse_param(stream, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(stream, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(stream, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(stream, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(stream, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(stream, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "checkerboard") {
        auto tvalue = pbrt_texture_checkerboard{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "dimension") {
                parse_param(stream, ptype, tvalue.dimension);
            } else if (pname == "tex1") {
                parse_param(stream, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(stream, ptype, tvalue.tex2);
            } else if (pname == "aamode") {
                parse_param(stream, ptype, tvalue.aamode);
            } else if (pname == "mapping") {
                parse_param(stream, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(stream, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(stream, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(stream, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(stream, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(stream, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "dots") {
        auto tvalue = pbrt_texture_dots{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "inside") {
                parse_param(stream, ptype, tvalue.inside);
            } else if (pname == "outside") {
                parse_param(stream, ptype, tvalue.outside);
            } else if (pname == "mapping") {
                parse_param(stream, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(stream, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(stream, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(stream, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(stream, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(stream, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "imagemap") {
        auto tvalue = pbrt_texture_imagemap{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "filename") {
                parse_param(stream, ptype, tvalue.filename);
            } else if (pname == "wrap") {
                parse_param(stream, ptype, tvalue.wrap);
            } else if (pname == "maxanisotropy") {
                parse_param(stream, ptype, tvalue.maxanisotropy);
            } else if (pname == "trilinear") {
                parse_param(stream, ptype, tvalue.trilinear);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "gamma") {
                parse_param(stream, ptype, tvalue.gamma);
            } else if (pname == "mapping") {
                parse_param(stream, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(stream, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(stream, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(stream, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(stream, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(stream, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_texture_mix{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "tex1") {
                parse_param(stream, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(stream, ptype, tvalue.tex2);
            } else if (pname == "amount") {
                parse_param(stream, ptype, tvalue.amount);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "scale") {
        auto tvalue = pbrt_texture_scale{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "tex1") {
                parse_param(stream, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(stream, ptype, tvalue.tex2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fbm") {
        auto tvalue = pbrt_texture_fbm{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "octaves") {
                parse_param(stream, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(stream, ptype, tvalue.roughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "wrinkled") {
        auto tvalue = pbrt_texture_wrinkled{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "octaves") {
                parse_param(stream, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(stream, ptype, tvalue.roughness);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "windy") {
        auto tvalue = pbrt_texture_windy{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "") {
                // TODO: missing params
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "marble") {
        auto tvalue = pbrt_texture_marble{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "octaves") {
                parse_param(stream, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(stream, ptype, tvalue.roughness);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "variation") {
                parse_param(stream, ptype, tvalue.variation);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uv") {
        auto tvalue = pbrt_texture_uv{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "mapping") {
                parse_param(stream, ptype, tvalue.mapping);
            } else if (pname == "uscale") {
                parse_param(stream, ptype, tvalue.uscale);
            } else if (pname == "vscale") {
                parse_param(stream, ptype, tvalue.vscale);
            } else if (pname == "udelta") {
                parse_param(stream, ptype, tvalue.udelta);
            } else if (pname == "vdelta") {
                parse_param(stream, ptype, tvalue.vdelta);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "v2") {
                parse_param(stream, ptype, tvalue.v2);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Texture " + type);
    }
}

// Get typename
static inline void parse_typeparam(pbrt_token_stream& stream, string& value) {
    save_stream_position(stream);
    value      = "";
    auto pname = ""s, ptype = ""s;
    while (is_param(stream) && value == "") {
        parse_nametype(stream, pname, ptype);
        if (pname == "type") {
            parse_param(stream, ptype, value);
        } else {
            skip_param(stream);
        }
    }
    if (value == "") throw pbrtio_error("type not found");
    restore_stream_position(stream);
}

// Parse Material
static inline void parse_pbrt_material(
    pbrt_token_stream& stream, const string& type, pbrt_material& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "matte") {
        auto tvalue = pbrt_material_matte{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "sigma") {
                parse_param(stream, ptype, tvalue.sigma);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mirror") {
        auto tvalue = pbrt_material_mirror{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kr") {
                parse_param(stream, ptype, tvalue.Kr);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "plastic") {
        auto tvalue = pbrt_material_plastic{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(stream, ptype, tvalue.Ks);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "metal") {
        auto tvalue = pbrt_material_metal{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "k") {
                parse_param(stream, ptype, tvalue.k);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "glass") {
        auto tvalue = pbrt_material_glass{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kr") {
                parse_param(stream, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(stream, ptype, tvalue.Kt);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "translucent") {
        auto tvalue = pbrt_material_translucent{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(stream, ptype, tvalue.Ks);
            } else if (pname == "reflect") {
                parse_param(stream, ptype, tvalue.reflect);
            } else if (pname == "transmit") {
                parse_param(stream, ptype, tvalue.transmit);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uber") {
        auto tvalue = pbrt_material_uber{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(stream, ptype, tvalue.Ks);
            } else if (pname == "Kr") {
                parse_param(stream, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(stream, ptype, tvalue.Kt);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "opacity") {
                parse_param(stream, ptype, tvalue.opacity);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disney") {
        auto tvalue = pbrt_material_disney{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "color") {
                parse_param(stream, ptype, tvalue.color);
            } else if (pname == "anisotropic") {
                parse_param(stream, ptype, tvalue.anisotropic);
            } else if (pname == "clearcoat") {
                parse_param(stream, ptype, tvalue.clearcoat);
            } else if (pname == "clearcoatgloss") {
                parse_param(stream, ptype, tvalue.clearcoatgloss);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "metallic") {
                parse_param(stream, ptype, tvalue.metallic);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "scatterdistance") {
                parse_param(stream, ptype, tvalue.scatterdistance);
            } else if (pname == "sheen") {
                parse_param(stream, ptype, tvalue.sheen);
            } else if (pname == "sheentint") {
                parse_param(stream, ptype, tvalue.sheentint);
            } else if (pname == "spectrans") {
                parse_param(stream, ptype, tvalue.spectrans);
            } else if (pname == "thin") {
                parse_param(stream, ptype, tvalue.thin);
            } else if (pname == "difftrans") {
                parse_param(stream, ptype, tvalue.difftrans);
            } else if (pname == "flatness") {
                parse_param(stream, ptype, tvalue.flatness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hair") {
        auto tvalue = pbrt_material_hair{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "color") {
                parse_param(stream, ptype, tvalue.color);
            } else if (pname == "sigma_a") {
                parse_param(stream, ptype, tvalue.sigma_a);
            } else if (pname == "eumelanin") {
                parse_param(stream, ptype, tvalue.eumelanin);
            } else if (pname == "pheomelanin") {
                parse_param(stream, ptype, tvalue.pheomelanin);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "beta_m") {
                parse_param(stream, ptype, tvalue.beta_m);
            } else if (pname == "beta_n") {
                parse_param(stream, ptype, tvalue.beta_n);
            } else if (pname == "alpha") {
                parse_param(stream, ptype, tvalue.alpha);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "kdsubsurface") {
        auto tvalue = pbrt_material_kdsubsurface{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "Kr") {
                parse_param(stream, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(stream, ptype, tvalue.Kt);
            } else if (pname == "mfp") {
                parse_param(stream, ptype, tvalue.mfp);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_material_mix{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "amount") {
                parse_param(stream, ptype, tvalue.amount);
            } else if (pname == "namedmaterial1") {
                parse_param(stream, ptype, tvalue.namedmaterial1);
            } else if (pname == "namedmaterial2") {
                parse_param(stream, ptype, tvalue.namedmaterial2);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fourier") {
        auto tvalue = pbrt_material_fourier{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "bsdffile") {
                parse_param(stream, ptype, tvalue.bsdffile);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "substrate") {
        auto tvalue = pbrt_material_substrate{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_param(stream, ptype, tvalue.Kd);
            } else if (pname == "Ks") {
                parse_param(stream, ptype, tvalue.Ks);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "subsurface") {
        auto tvalue = pbrt_material_subsurface{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "name") {
                parse_param(stream, ptype, tvalue.name);
            } else if (pname == "sigma_a") {
                parse_param(stream, ptype, tvalue.sigma_a);
            } else if (pname == "sigma_prime_s") {
                parse_param(stream, ptype, tvalue.sigma_prime_s);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "eta") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "index") {
                parse_param(stream, ptype, tvalue.eta);
            } else if (pname == "Kr") {
                parse_param(stream, ptype, tvalue.Kr);
            } else if (pname == "Kt") {
                parse_param(stream, ptype, tvalue.Kt);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_param(stream, ptype, tvalue.uroughness);
            } else if (pname == "vroughness") {
                parse_param(stream, ptype, tvalue.vroughness);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_param(stream, ptype, tvalue.bumpmap);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Material " + type);
    }
}

// Parse Shape
static inline void parse_pbrt_shape(
    pbrt_token_stream& stream, const string& type, pbrt_shape& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "trianglemesh") {
        auto tvalue = pbrt_shape_trianglemesh{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "indices") {
                parse_param(stream, ptype, tvalue.indices);
            } else if (pname == "P") {
                parse_param(stream, ptype, tvalue.P);
            } else if (pname == "N") {
                parse_param(stream, ptype, tvalue.N);
            } else if (pname == "S") {
                parse_param(stream, ptype, tvalue.S);
            } else if (pname == "uv") {
                parse_param(stream, ptype, tvalue.uv);
            } else if (pname == "st") {
                parse_param(stream, ptype, tvalue.uv);
            } else if (pname == "alpha") {
                parse_param(stream, ptype, tvalue.alpha);
            } else if (pname == "shadowalpha") {
                parse_param(stream, ptype, tvalue.shadowalpha);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "plymesh") {
        auto tvalue = pbrt_shape_plymesh{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "filename") {
                parse_param(stream, ptype, tvalue.filename);
            } else if (pname == "alpha") {
                parse_param(stream, ptype, tvalue.alpha);
            } else if (pname == "shadowalpha") {
                parse_param(stream, ptype, tvalue.shadowalpha);
            } else if (pname == "discarddegenerateUVs") {
                // hack for some files
                auto value = false;
                parse_param(stream, ptype, value);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "curve") {
        auto tvalue = pbrt_shape_curve{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "P") {
                parse_param(stream, ptype, tvalue.P);
            } else if (pname == "N") {
                parse_param(stream, ptype, tvalue.N);
            } else if (pname == "basis") {
                parse_param(stream, ptype, tvalue.basis);
            } else if (pname == "degree") {
                parse_param(stream, ptype, tvalue.degree);
            } else if (pname == "type") {
                parse_param(stream, ptype, tvalue.type);
            } else if (pname == "width") {
                auto width = 1.0f;
                parse_param(stream, ptype, width);
                tvalue.width0 = width;
                tvalue.width1 = width;
            } else if (pname == "width0") {
                parse_param(stream, ptype, tvalue.width0);
            } else if (pname == "width1") {
                parse_param(stream, ptype, tvalue.width1);
            } else if (pname == "splitdepth") {
                parse_param(stream, ptype, tvalue.splitdepth);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "loopsubdiv") {
        auto tvalue = pbrt_shape_loopsubdiv{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "indices") {
                parse_param(stream, ptype, tvalue.indices);
            } else if (pname == "P") {
                parse_param(stream, ptype, tvalue.P);
            } else if (pname == "levels") {
                parse_param(stream, ptype, tvalue.levels);
            } else if (pname == "nlevels") {
                parse_param(stream, ptype, tvalue.levels);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "nurbs") {
        auto tvalue = pbrt_shape_nurbs{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "nu") {
                parse_param(stream, ptype, tvalue.nu);
            } else if (pname == "nv") {
                parse_param(stream, ptype, tvalue.nv);
            } else if (pname == "uknots") {
                parse_param(stream, ptype, tvalue.uknots);
            } else if (pname == "vknots") {
                parse_param(stream, ptype, tvalue.vknots);
            } else if (pname == "u0") {
                parse_param(stream, ptype, tvalue.u0);
            } else if (pname == "v0") {
                parse_param(stream, ptype, tvalue.v0);
            } else if (pname == "u1") {
                parse_param(stream, ptype, tvalue.u1);
            } else if (pname == "v1") {
                parse_param(stream, ptype, tvalue.v1);
            } else if (pname == "P") {
                parse_param(stream, ptype, tvalue.P);
            } else if (pname == "Pw") {
                parse_param(stream, ptype, tvalue.Pw);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sphere") {
        auto tvalue = pbrt_shape_sphere{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(stream, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(stream, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disk") {
        auto tvalue = pbrt_shape_disk{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "height") {
                parse_param(stream, ptype, tvalue.height);
            } else if (pname == "innerradius") {
                parse_param(stream, ptype, tvalue.innerradius);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cone") {
        auto tvalue = pbrt_shape_cone{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "height") {
                parse_param(stream, ptype, tvalue.height);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cylinder") {
        auto tvalue = pbrt_shape_cylinder{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(stream, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(stream, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hyperboloid") {
        auto tvalue = pbrt_shape_hyperboloid{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "p1") {
                parse_param(stream, ptype, tvalue.p1);
            } else if (pname == "p2") {
                parse_param(stream, ptype, tvalue.p2);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "paraboloid") {
        auto tvalue = pbrt_shape_paraboloid{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "zmin") {
                parse_param(stream, ptype, tvalue.zmin);
            } else if (pname == "zmax") {
                parse_param(stream, ptype, tvalue.zmax);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heightfield") {
        auto tvalue = pbrt_shape_heightfield{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "nu") {
                parse_param(stream, ptype, tvalue.nu);
            } else if (pname == "nv") {
                parse_param(stream, ptype, tvalue.nv);
            } else if (pname == "Pz") {
                parse_param(stream, ptype, tvalue.Pz);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Shape " + type);
    }
}

// Parse AreaLightSource
static inline void parse_pbrt_arealight(
    pbrt_token_stream& stream, const string& type, pbrt_arealight& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "diffuse") {
        auto tvalue = pbrt_arealight_diffuse{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "L") {
                parse_param(stream, ptype, tvalue.L);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "twosided") {
                parse_param(stream, ptype, tvalue.twosided);
            } else if (pname == "samples") {
                parse_param(stream, ptype, tvalue.samples);
            } else if (pname == "nsamples") {
                parse_param(stream, ptype, tvalue.samples);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse LightSource
static inline void parse_pbrt_light(
    pbrt_token_stream& stream, const string& type, pbrt_light& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "distant") {
        auto tvalue = pbrt_light_distant{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(stream, ptype, tvalue.L);
            } else if (pname == "from") {
                parse_param(stream, ptype, tvalue.from);
            } else if (pname == "to") {
                parse_param(stream, ptype, tvalue.to);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "goniometric") {
        auto tvalue = pbrt_light_goniometric{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(stream, ptype, tvalue.I);
            } else if (pname == "mapname") {
                parse_param(stream, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "infinite") {
        auto tvalue = pbrt_light_infinite{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(stream, ptype, tvalue.L);
            } else if (pname == "samples") {
                parse_param(stream, ptype, tvalue.samples);
            } else if (pname == "nsamples") {
                parse_param(stream, ptype, tvalue.samples);
            } else if (pname == "mapname") {
                parse_param(stream, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "distant") {
        auto tvalue = pbrt_light_distant{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(stream, ptype, tvalue.L);
            } else if (pname == "from") {
                parse_param(stream, ptype, tvalue.from);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "projection") {
        auto tvalue = pbrt_light_projection{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(stream, ptype, tvalue.I);
            } else if (pname == "fov") {
                parse_param(stream, ptype, tvalue.fov);
            } else if (pname == "mapname") {
                parse_param(stream, ptype, tvalue.mapname);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "spot") {
        auto tvalue = pbrt_light_spot{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(stream, ptype, tvalue.I);
            } else if (pname == "from") {
                parse_param(stream, ptype, tvalue.from);
            } else if (pname == "to") {
                parse_param(stream, ptype, tvalue.to);
            } else if (pname == "coneangle") {
                parse_param(stream, ptype, tvalue.coneangle);
            } else if (pname == "conedeltaangle") {
                parse_param(stream, ptype, tvalue.conedeltaangle);
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Film " + type);
    }
}

// Parse Medium
static inline void parse_pbrt_medium(
    pbrt_token_stream& stream, const string& type, pbrt_medium& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "homogeneous") {
        auto tvalue = pbrt_medium_homogeneous{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "sigma_a") {
                parse_param(stream, ptype, tvalue.sigma_a);
            } else if (pname == "sigma_s") {
                parse_param(stream, ptype, tvalue.sigma_s);
            } else if (pname == "preset") {
                parse_param(stream, ptype, tvalue.preset);
            } else if (pname == "g") {
                parse_param(stream, ptype, tvalue.g);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heterogeneous") {
        auto tvalue = pbrt_medium_heterogeneous{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "sigma_a") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "sigma_s") {
                parse_param(stream, ptype, tvalue.sigma_s);
            } else if (pname == "preset") {
                parse_param(stream, ptype, tvalue.preset);
            } else if (pname == "g") {
                parse_param(stream, ptype, tvalue.g);
            } else if (pname == "p0") {
                parse_param(stream, ptype, tvalue.p0);
            } else if (pname == "p1") {
                parse_param(stream, ptype, tvalue.p1);
            } else if (pname == "nx") {
                parse_param(stream, ptype, tvalue.nx);
            } else if (pname == "ny") {
                parse_param(stream, ptype, tvalue.ny);
            } else if (pname == "nz") {
                parse_param(stream, ptype, tvalue.nz);
            } else if (pname == "density") {
                parse_param(stream, ptype, tvalue.density);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw pbrtio_error("inconsistent types");
            } else {
                throw pbrtio_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw pbrtio_error("unknown Medium " + type);
    }
}

// Load pbrt scene
template <typename Callbacks>
inline void load_pbrt(
    const string& filename, Callbacks& cb, const load_pbrt_options& options) {
    // start laoding files
    auto streams = vector<pbrt_token_stream>{};
    init_token_stream(streams);
    load_token_stream(filename, streams);

    // parsing stack
    auto stack    = vector<pbrt_context>{{}};
    auto object   = pbrt_object{};
    auto coordsys = unordered_map<string, affine3f>{};

    // parse command by command
    auto cmd = ""s;
    while (!streams.empty() && !is_empty(streams.back())) {
        // get command
        auto& stream = streams.back();
        parse_command(stream, cmd);
        if (cmd == "WorldBegin") {
            stack.push_back({});
        } else if (cmd == "WorldEnd") {
            stack.pop_back();
            if (stack.size() != 1) throw pbrtio_error("bad stack");
        } else if (cmd == "AttributeBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "AttributeEnd") {
            stack.pop_back();
        } else if (cmd == "TransformBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "TransformEnd") {
            stack.pop_back();
        } else if (cmd == "ObjectBegin") {
            parse_value(stream, object.name);
            cb.begin_object(object, stack.back());
        } else if (cmd == "ObjectEnd") {
            cb.end_object(object, stack.back());
            object = {};
        } else if (cmd == "ObjectInstance") {
            auto value = pbrt_object{};
            parse_value(stream, value.name);
            cb.object_instance(value, stack.back());
        } else if (cmd == "Transform") {
            auto xf = identity_mat4f;
            parse_param(stream, xf);
            stack.back().frame = (affine3f)xf;
        } else if (cmd == "ConcatTransform") {
            auto xf = identity_mat4f;
            parse_param(stream, xf);
            stack.back().frame = stack.back().frame * (affine3f)xf;
        } else if (cmd == "Scale") {
            auto v = zero3f;
            parse_param(stream, v);
            stack.back().frame = stack.back().frame *
                                 (const affine3f&)make_scaling_frame(v);
        } else if (cmd == "Translate") {
            auto v = zero3f;
            parse_param(stream, v);
            stack.back().frame = stack.back().frame *
                                 (const affine3f&)make_translation_frame(v);
        } else if (cmd == "Rotate") {
            auto v = zero4f;
            parse_param(stream, v);
            stack.back().frame = stack.back().frame *
                                 (const affine3f&)make_rotation_frame(
                                     vec3f{v.y, v.z, v.w}, radians(v.x));
        } else if (cmd == "LookAt") {
            auto from = zero3f, to = zero3f, up = zero3f;
            parse_param(stream, from);
            parse_param(stream, to);
            parse_param(stream, up);
            // from pbrt parser
            auto frame = make_lookat_frame(from, to, up, true);
            // frame.z = normalize(to-from);
            // frame.x = normalize(cross(frame.z,up));
            // frame.y = cross(frame.x,frame.z);
            // frame.o    = from;
            stack.back().frame = stack.back().frame *
                                 (const affine3f&)inverse(frame);
            // stack.back().focus = length(m.x - m.y);
        } else if (cmd == "ReverseOrientation") {
            stack.back().reverse = !stack.back().reverse;
        } else if (cmd == "CoordinateSystem") {
            auto name = ""s;
            parse_value(stream, name);
            coordsys[name] = stack.back().frame;
        } else if (cmd == "CoordSysTransform") {
            auto name = ""s;
            parse_value(stream, name);
            if (coordsys.find(name) != coordsys.end()) {
                stack.back().frame = coordsys.at(name);
            }
        } else if (cmd == "Integrator") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_integrator{};
            parse_pbrt_integrator(stream, type, value);
            cb.integrator(value, stack.back());
        } else if (cmd == "Sampler") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_sampler{};
            parse_pbrt_sampler(stream, type, value);
            cb.sampler(value, stack.back());
        } else if (cmd == "PixelFilter") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_filter{};
            parse_pbrt_filter(stream, type, value);
            cb.filter(value, stack.back());
        } else if (cmd == "Film") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_film{};
            parse_pbrt_film(stream, type, value);
            cb.film(value, stack.back());
        } else if (cmd == "Camera") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_camera{};
            parse_pbrt_camera(stream, type, value);
            cb.camera(value, stack.back());
        } else if (cmd == "Texture") {
            auto name = ""s, comptype = ""s, type = ""s;
            parse_value(stream, name);
            parse_value(stream, comptype);
            parse_value(stream, type);
            auto value = pbrt_texture{};
            parse_pbrt_texture(stream, type, value);
            cb.texture(value, name, stack.back());
        } else if (cmd == "Material") {
            static auto material_id = 0;
            auto        type        = ""s;
            parse_value(stream, type);
            if (type == "") {
                stack.back().material = "";
            } else {
                auto value = pbrt_material{};
                auto name = "unnamed_material_" + std::to_string(material_id++);
                parse_pbrt_material(stream, type, value);
                stack.back().material = name;
                cb.material(value, name, stack.back());
            }
        } else if (cmd == "MakeNamedMaterial") {
            auto name = ""s, type = ""s;
            parse_value(stream, name);
            parse_typeparam(stream, type);
            auto value = pbrt_material{};
            parse_pbrt_material(stream, type, value);
            cb.material(value, name, stack.back());
        } else if (cmd == "NamedMaterial") {
            auto name = ""s;
            parse_value(stream, name);
            stack.back().material = name;
        } else if (cmd == "Shape") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_shape{};
            parse_pbrt_shape(stream, type, value);
            cb.shape(value, stack.back());
        } else if (cmd == "AreaLightSource") {
            auto type = ""s;
            parse_value(stream, type);
            static auto material_id = 0;
            auto name  = "unnamed_arealight_" + std::to_string(material_id++);
            auto value = pbrt_arealight{};
            parse_pbrt_arealight(stream, type, value);
            stack.back().arealight = name;
            cb.arealight(value, name, stack.back());
        } else if (cmd == "LightSource") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_light{};
            parse_pbrt_light(stream, type, value);
            cb.light(value, stack.back());
        } else if (cmd == "MakeNamedMedium") {
            auto name = ""s, type = ""s;
            parse_value(stream, name);
            parse_typeparam(stream, type);
            auto value = pbrt_medium{};
            parse_pbrt_medium(stream, type, value);
            cb.medium(value, name, stack.back());
        } else if (cmd == "MediumInterface") {
            auto interior = ""s, exterior = ""s;
            parse_value(stream, interior);
            parse_value(stream, exterior);
            stack.back().medium_interior = interior;
            stack.back().medium_exterior = exterior;
        } else if (cmd == "Include") {
            auto inputname = ""s;
            parse_value(stream, inputname);
            load_token_stream(get_dirname(filename) + inputname, streams);
        } else {
            throw pbrtio_error("unknown command " + cmd);
        }
        // go to next file if needed
        skip_whitespace_or_comment_to_next_file(streams);
    }
}

}  // namespace yocto

#endif
