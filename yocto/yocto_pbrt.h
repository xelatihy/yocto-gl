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
template <typename T, int N>
struct pbrt_spectrum;

// pbrt pbrt_spectrum as rgb color
template <typename T>
struct pbrt_spectrum<T, 3> {
    union {
        struct {
            T x, y, z;
        };
        T elems[3];
    };

    constexpr pbrt_spectrum() : x{0}, y{0}, z{0} {}
    constexpr pbrt_spectrum(T x, T y, T z) : x{x}, y{y}, z{z} {}
    constexpr explicit pbrt_spectrum(T v) : x{v}, y{v}, z{v} {}
    constexpr explicit operator vec<T, 3>() const { return {x, y, z}; };

    constexpr T&       operator[](int i) { return elems[i]; }
    constexpr const T& operator[](int i) const { return elems[i]; }
};

// typedefs
using pbrt_spectrum3f = pbrt_spectrum<float, 3>;

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
using pbrt_camera = variant<pbrt_perspective_camera, pbrt_orthographic_camera,
    pbrt_environment_camera, pbrt_realistic_camera>;

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
using pbrt_sampler = variant<pbrt_random_sampler, pbrt_halton_sampler,
    pbrt_sobol_sampler, pbrt_zerotwosequence_sampler, pbrt_maxmindist_sampler,
    pbrt_stratified_sampler>;

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
using pbrt_filter = variant<pbrt_box_filter, pbrt_gaussian_filter,
    pbrt_mitchell_filter, pbrt_sinc_filter, pbrt_triangle_filter>;

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
using pbrt_integrator = variant<pbrt_path_integrator, pbrt_volpath_integrator,
    pbrt_bdpt_integrator, pbrt_directlighting_integrator, pbrt_mlt_integrator,
    pbrt_sppm_integrator, pbrt_whitted_integrator>;

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
using pbrt_accelerator = variant<pbrt_bvh_accelerator, pbrt_kdtree_accelerator>;

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
    pbrt_textured<pbrt_spectrum3f> Kd      = {0.5f, 0.5f, 0.5f};
    pbrt_textured<float>           sigma   = 0;
    pbrt_textured<float>           bumpmap = 0;
};
struct pbrt_mirror_material {
    pbrt_textured<pbrt_spectrum3f> Kr      = {0.9f, 0.9f, 0.9f};
    pbrt_textured<float>           bumpmap = 0;
};
struct pbrt_plastic_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0.25f, 0.25f, 0.25f};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0.25f, 0.25f, 0.25f};
    pbrt_textured<float>           uroughness     = 0.1f;
    pbrt_textured<float>           vroughness     = 0.1f;
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
    pbrt_textured<float>           eta            = 1;
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_translucent_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> reflect        = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> transmit       = {0, 0, 0};
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_uber_material {
    pbrt_textured<pbrt_spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Kr             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Kt             = {0, 0, 0};
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
    pbrt_textured<float>           eta            = 1;
    pbrt_textured<pbrt_spectrum3f> opacity        = {1, 1, 1};
    bool                           remaproughness = true;
    pbrt_textured<float>           bumpmap        = 0;
};
struct pbrt_disney_material {
    pbrt_textured<pbrt_spectrum3f> color           = {0.5f, 0.5f, 0.5f};
    pbrt_textured<float>           anisotropic     = 0;
    pbrt_textured<float>           clearcoat       = 0;
    pbrt_textured<float>           clearcoatgloss  = 1;
    pbrt_textured<float>           eta             = 1.5f;
    pbrt_textured<float>           metallic        = 0;
    pbrt_textured<float>           uroughness      = 0;
    pbrt_textured<float>           vroughness      = 0;
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
    pbrt_textured<pbrt_spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> mfp            = {1, 1, 1};
    pbrt_textured<float>           eta            = 1;
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
    pbrt_textured<pbrt_spectrum3f> Kd             = {0, 0, 0};
    pbrt_textured<pbrt_spectrum3f> Ks             = {0, 0, 0};
    pbrt_textured<float>           uroughness     = 0;
    pbrt_textured<float>           vroughness     = 0;
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
using pbrt_arealight = variant<pbrt_none_arealight, pbrt_diffuse_arealight>;

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
using pbrt_medium = variant<pbrt_homogeneous_medium, pbrt_heterogeneous_medium>;

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
    void sampler(const pbrt_sampler& value, const pbrt_context& ctx) {}
    void integrator(const pbrt_integrator& value, const pbrt_context& ctx) {}
    void accelerator(const pbrt_accelerator& value, const pbrt_context& ctx) {}
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

// Load pbrt params
struct pbrt_params {
    bool geometry_only = false;
    bool flip_texcoord = true;
};

// Load pbrt scene
template <typename Callbacks>
inline void load_pbrt(
    const string& filename, Callbacks& cb, const pbrt_params& params = {});

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
struct pbrt_stream {
    string      buffer;
    string_view str;
    string_view saved;
};

// skip white space or comment
static inline void skip_whitespace_or_comment(pbrt_stream& stream) {
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
static inline void parse_value(pbrt_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.front() != '"') {
        throw io_error("bad string");
    }
    str.remove_prefix(1);
    auto pos = str.find('"');
    if (pos == string_view::npos) {
        throw io_error("bad string");
    }
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
}

// parse a quoted string
static inline void parse_command(pbrt_stream& stream, string& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (!std::isalpha((int)str.front())) {
        throw io_error("bad command");
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
static inline void parse_value(pbrt_stream& stream, float& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw io_error("number expected");
    auto next = (char*)nullptr;
    value     = strtof(str.data(), &next);
    if (str.data() == next) throw io_error("number expected");
    str.remove_prefix(next - str.data());
}

// parse a number
static inline void parse_value(pbrt_stream& stream, int& value) {
    skip_whitespace_or_comment(stream);
    auto& str = stream.str;
    if (str.empty()) throw io_error("number expected");
    auto next = (char*)nullptr;
    value     = strtol(str.data(), &next, 10);
    if (str.data() == next) throw io_error("number expected");
    str.remove_prefix(next - str.data());
}
static inline void parse_value(pbrt_stream& stream, bool& value) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    if (value_name == "true") {
        value = true;
    } else if (value_name == "false") {
        value = false;
    } else {
        throw io_error("expected boolean");
    }
}
template <typename T>
static inline void parse_value(
    pbrt_stream& stream, T& value, unordered_map<string, T>& value_names) {
    auto value_name = ""s;
    parse_value(stream, value_name);
    try {
        value = value_names.at(value_name);
    } catch (std::out_of_range&) {
        throw io_error("expected enum value");
    }
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_bilerp_texture::mapping_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_bilerp_texture::mapping_type>{
            {"uv", pbrt_bilerp_texture::mapping_type::uv},
            {"spherical", pbrt_bilerp_texture::mapping_type::spherical},
            {"cylindrical", pbrt_bilerp_texture::mapping_type::cylindrical},
            {"planar", pbrt_bilerp_texture::mapping_type::planar},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_checkerboard_texture::mapping_type& value) {
    return parse_value(stream, (pbrt_bilerp_texture::mapping_type&)value);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_dots_texture::mapping_type& value) {
    return parse_value(stream, (pbrt_bilerp_texture::mapping_type&)value);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_imagemap_texture::mapping_type& value) {
    return parse_value(stream, (pbrt_bilerp_texture::mapping_type&)value);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_uv_texture::mapping_type& value) {
    return parse_value(stream, (pbrt_bilerp_texture::mapping_type&)value);
}

static inline void parse_value(
    pbrt_stream& stream, pbrt_checkerboard_texture::aamode_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_checkerboard_texture::aamode_type>{
            {"closedform", pbrt_checkerboard_texture::aamode_type::closedform},
            {"none", pbrt_checkerboard_texture::aamode_type::none},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_imagemap_texture::wrap_type& value) {
    static auto value_names =
        unordered_map<string, pbrt_imagemap_texture::wrap_type>{
            {"repeat", pbrt_imagemap_texture::wrap_type::repeat},
            {"clamp", pbrt_imagemap_texture::wrap_type::clamp},
            {"black", pbrt_imagemap_texture::wrap_type::black},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_curve_shape::basis_t& value) {
    static auto value_names = unordered_map<string, pbrt_curve_shape::basis_t>{
        {"bezier", pbrt_curve_shape::basis_t::bezier},
        {"bspline", pbrt_curve_shape::basis_t::bspline},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_curve_shape::type_t& value) {
    static auto value_names = unordered_map<string, pbrt_curve_shape::type_t>{
        {"flat", pbrt_curve_shape::type_t::flat},
        {"cylinder", pbrt_curve_shape::type_t::cylinder},
        {"ribbon", pbrt_curve_shape::type_t::ribbon},
    };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_bvh_accelerator::splitmethod_t& value) {
    static auto value_names =
        unordered_map<string, pbrt_bvh_accelerator::splitmethod_t>{
            {"sah", pbrt_bvh_accelerator::splitmethod_t::sah},
            {"equal", pbrt_bvh_accelerator::splitmethod_t::equal},
            {"middle", pbrt_bvh_accelerator::splitmethod_t::middle},
            {"hlbvh", pbrt_bvh_accelerator::splitmethod_t::hlbvh},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_path_integrator::lightsamplestrategy_t& value) {
    static auto value_names =
        unordered_map<string, pbrt_path_integrator::lightsamplestrategy_t>{
            {"power", pbrt_path_integrator::lightsamplestrategy_t::power},
            {"spatial", pbrt_path_integrator::lightsamplestrategy_t::spatial},
            {"uniform", pbrt_path_integrator::lightsamplestrategy_t::uniform},
        };
    return parse_value(stream, value, value_names);
}
static inline void parse_value(pbrt_stream&         stream,
    pbrt_volpath_integrator::lightsamplestrategy_t& value) {
    return parse_value(
        stream, (pbrt_path_integrator::lightsamplestrategy_t&)value);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_bdpt_integrator::lightsamplestrategy_t& value) {
    return parse_value(
        stream, (pbrt_path_integrator::lightsamplestrategy_t&)value);
}
static inline void parse_value(
    pbrt_stream& stream, pbrt_directlighting_integrator::strategy_t& value) {
    static auto value_names =
        unordered_map<string, pbrt_directlighting_integrator::strategy_t>{
            {"all", pbrt_directlighting_integrator::strategy_t::all},
            {"one", pbrt_directlighting_integrator::strategy_t::one},
        };
    return parse_value(stream, value, value_names);
}

// parse a vec type
template <typename T, int N>
static inline void parse_value(pbrt_stream& stream, vec<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}
template <typename T, int N, int M>
static inline void parse_value(pbrt_stream& stream, mat<T, N, M>& value) {
    for (auto i = 0; i < M; i++) parse_value(stream, value[i]);
}
template <typename T>
static inline void parse_value(pbrt_stream& stream, bbox<T, 2>& value) {
    parse_value(stream, value[0][0]);
    parse_value(stream, value[1][0]);
    parse_value(stream, value[0][1]);
    parse_value(stream, value[1][1]);
}
template <typename T, int N>
static inline void parse_value(
    pbrt_stream& stream, pbrt_spectrum<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(stream, value[i]);
}

// Check next
static inline bool is_empty(pbrt_stream& stream) {
    skip_whitespace_or_comment(stream);
    return stream.str.empty();
}
static inline bool is_string(pbrt_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '"';
}
static inline bool is_open_bracket(pbrt_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == '[';
}
static inline bool is_close_bracket(pbrt_stream& stream) {
    skip_whitespace_or_comment(stream);
    return !stream.str.empty() && stream.str.front() == ']';
}
static inline bool is_param(pbrt_stream& stream) {
    skip_whitespace_or_comment(stream);
    return is_string(stream);
}

// parse a quoted string
static inline void parse_nametype(
    pbrt_stream& stream, string& name, string& type) {
    auto value = ""s;
    parse_value(stream, value);
    auto str  = string_view{value};
    auto pos1 = str.find(' ');
    if (pos1 == string_view::npos) {
        throw io_error("bad type " + value);
    }
    type = string(str.substr(0, pos1));
    str.remove_prefix(pos1);
    auto pos2 = str.find_first_not_of(' ');
    if (pos2 == string_view::npos) {
        throw io_error("bad type " + value);
    }
    str.remove_prefix(pos2);
    name = string(str);
}

static inline void skip_open_bracket(pbrt_stream& stream) {
    if (!is_open_bracket(stream)) throw io_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}
static inline void skip_close_bracket(pbrt_stream& stream) {
    if (!is_close_bracket(stream)) throw io_error("expected bracket");
    stream.str.remove_prefix(1);
    skip_whitespace_or_comment(stream);
}

template <typename T>
static inline void parse_param(pbrt_stream& stream, T& value) {
    auto has_brackets = is_open_bracket(stream);
    if (has_brackets) skip_open_bracket(stream);
    parse_value(stream, value);
    if (has_brackets) skip_close_bracket(stream);
}

template <typename T>
static inline void parse_param(pbrt_stream& stream, vector<T>& values) {
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
    } else if constexpr (std::is_same<T, pbrt_spectrum3f>::value) {
        return type == "rgb" || type == "pbrt_spectrum" || type == "blackbody";
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
    pbrt_stream& stream, const string& type, T& value) {
    if (!is_type_compatible<T>(type)) {
        throw io_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

static inline pair<vec3f, vec3f> pbrt_get_element_etak(const string& name) {
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

template <typename T>
static inline void parse_param(
    pbrt_stream& stream, const string& type, pbrt_spectrum<T, 3>& value) {
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
    } else if (type == "pbrt_spectrum" && is_string(stream)) {
        if (verbose) printf("pbrt_spectrum  not well supported\n");
        auto filename = ""s;
        parse_param(stream, filename);
        filename = get_filename(filename);
        if (get_extension(filename) == "spd") {
            filename = filename.substr(0, filename.size() - 4);
            if (filename == "SHPS") {
                value = {1, 1, 1};
            } else if (get_extension(filename) == "eta") {
                auto eta = pbrt_get_element_etak(
                    filename.substr(0, filename.length() - 4))
                               .first;
                value = {eta.x, eta.y, eta.z};
            } else if (get_extension(filename) == "k") {
                auto k = pbrt_get_element_etak(
                    filename.substr(0, filename.length() - 2))
                             .second;
                value = {k.x, k.y, k.z};
            } else {
                throw io_error("unknown pbrt_spectrum file " + filename);
            }
        } else {
            throw io_error("unsupported pbrt_spectrum format");
            // value = {1, 0, 0};
        }
    } else if (type == "pbrt_spectrum" && !is_string(stream)) {
        if (verbose) printf("pbrt_spectrum  not well supported\n");
        auto values = vector<float>{};
        parse_param(stream, values);
        value = {1, 0, 0};
    } else {
        throw io_error("unsupported pbrt_spectrum type");
    }
}

template <typename T>
static inline void parse_param(
    pbrt_stream& stream, const string& type, vector<T>& value) {
    if (!is_type_compatible<T>(type)) {
        throw io_error("incompatible type " + type);
    }
    parse_param(stream, value);
}

template <typename T>
static inline void parse_param(
    pbrt_stream& stream, const string& type, pbrt_textured<T>& value) {
    if (type == "texture") {
        parse_param(stream, value.texture);
    } else {
        parse_param(stream, type, value.value);
    }
}

static inline void skip_value(pbrt_stream& stream) {
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

static inline void skip_param(pbrt_stream& stream) {
    if (is_open_bracket(stream)) {
        skip_open_bracket(stream);
        while (!is_close_bracket(stream)) skip_value(stream);
        skip_close_bracket(stream);
    } else {
        skip_value(stream);
    }
}

static inline void save_stream_position(pbrt_stream& stream) {
    stream.saved = stream.str;
}
static inline void restore_stream_position(pbrt_stream& stream) {
    stream.str = stream.saved;
}

// Load a token stream
static inline void load_stream(const string& filename, pbrt_stream& stream) {
    load_text(filename, stream.buffer);
    stream.str = stream.buffer;
}

// operations on token stacks
static inline void init_stream(vector<pbrt_stream>& stream) {
    stream.reserve(100);
}
static inline void load_stream(
    const string& filename, vector<pbrt_stream>& stream) {
    stream.emplace_back();
    load_stream(filename, stream.back());
}

// Skip whitespace
static inline void skip_whitespace_or_comment_to_next_file(
    vector<pbrt_stream>& stream) {
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

// Parse Accelerator
static inline void parse_accelerator(
    pbrt_stream& stream, const string& type, pbrt_accelerator& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "bvh") {
        auto tvalue = pbrt_bvh_accelerator{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxnodeprims") {
                parse_param(stream, ptype, tvalue.maxnodeprims);
            } else if (pname == "splitmethod") {
                parse_param(stream, ptype, tvalue.splitmethod);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "kdtree") {
        auto tvalue = pbrt_kdtree_accelerator{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "intersectcost") {
                parse_param(stream, ptype, tvalue.intersectcost);
            } else if (pname == "traversalcost") {
                parse_param(stream, ptype, tvalue.traversalcost);
            } else if (pname == "emptybonus") {
                parse_param(stream, ptype, tvalue.emptybonus);
            } else if (pname == "maxprims") {
                parse_param(stream, ptype, tvalue.maxprims);
            } else if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Accelerator " + type);
    }
}

// Parse Integrator
static inline void parse_integrator(
    pbrt_stream& stream, const string& type, pbrt_integrator& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "path") {
        auto tvalue = pbrt_path_integrator{};
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
                throw io_error("unknown parameter " + pname);
            }
            // parse_optional_param(stream, "lightsamplestrategy",
            // tvalue.lightsamplestrategy); // TODO: enums
        }
        value = tvalue;
    } else if (type == "volpath") {
        auto tvalue = pbrt_volpath_integrator{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "directlighting") {
        auto tvalue = pbrt_directlighting_integrator{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "strategy") {
                parse_param(stream, ptype, tvalue.strategy);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bdpt") {
        auto tvalue = pbrt_bdpt_integrator{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mlt") {
        auto tvalue = pbrt_mlt_integrator{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sppm") {
        auto tvalue = pbrt_sppm_integrator{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else if (pname == "iterations") {
                parse_param(stream, ptype, tvalue.iterations);
            } else if (pname == "numiterations") {
                parse_param(stream, ptype, tvalue.iterations);
            } else if (pname == "photonsperiteration") {
                parse_param(stream, ptype, tvalue.photonsperiteration);
            } else if (pname == "imagewritefrequency") {
                parse_param(stream, ptype, tvalue.imagewritefrequency);
            } else if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "whitted") {
        auto tvalue = pbrt_whitted_integrator{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "maxdepth") {
                parse_param(stream, ptype, tvalue.maxdepth);
            } else if (pname == "pixelbounds") {
                parse_param(stream, ptype, tvalue.pixelbounds);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Integrator " + type);
    }
}

// Parse Sampler
static inline void parse_sampler(
    pbrt_stream& stream, const string& type, pbrt_sampler& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "random") {
        auto tvalue = pbrt_random_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "halton") {
        auto tvalue = pbrt_halton_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sobol") {
        auto tvalue = pbrt_sobol_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "02sequence") {
        auto tvalue = pbrt_zerotwosequence_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "lowdiscrepancy") {
        auto tvalue = pbrt_zerotwosequence_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "maxmindist") {
        auto tvalue = pbrt_maxmindist_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "pixelsamples") {
                parse_param(stream, ptype, tvalue.pixelsamples);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "stratified") {
        auto tvalue = pbrt_stratified_sampler{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xsamples") {
                parse_param(stream, ptype, tvalue.xsamples);
            } else if (pname == "ysamples") {
                parse_param(stream, ptype, tvalue.ysamples);
            } else if (pname == "jitter") {
                parse_param(stream, ptype, tvalue.jitter);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Sampler " + type);
    }
}

// Parse Filter
static inline void parse_filter(
    pbrt_stream& stream, const string& type, pbrt_filter& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "box") {
        auto tvalue = pbrt_box_filter{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "gaussian") {
        auto tvalue = pbrt_gaussian_filter{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else if (pname == "alpha") {
                parse_param(stream, ptype, tvalue.alpha);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mitchell") {
        auto tvalue = pbrt_mitchell_filter{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sinc") {
        auto tvalue = pbrt_sinc_filter{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else if (pname == "tau") {
                parse_param(stream, ptype, tvalue.tau);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "triangle") {
        auto tvalue = pbrt_triangle_filter{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "xwidth") {
                parse_param(stream, ptype, tvalue.xwidth);
            } else if (pname == "ywidth") {
                parse_param(stream, ptype, tvalue.ywidth);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown PixelFilter " + type);
    }
}

// Parse Filter
static inline void parse_film(
    pbrt_stream& stream, const string& type, pbrt_film& value) {
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Film " + type);
    }
}

// Parse Camera
static inline void parse_camera(
    pbrt_stream& stream, const string& type, pbrt_camera& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "perspective") {
        auto tvalue = pbrt_perspective_camera{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "orthographic") {
        auto tvalue = pbrt_orthographic_camera{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "environment") {
        auto tvalue = pbrt_environment_camera{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "shutteropen") {
                parse_param(stream, ptype, tvalue.shutteropen);
            } else if (pname == "shutterclose") {
                parse_param(stream, ptype, tvalue.shutterclose);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "realistic") {
        auto tvalue = pbrt_realistic_camera{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "lensfile") {
                parse_param(stream, ptype, tvalue.lensfile);
                // example: wide.22mm.dat
                auto lensfile = get_filename(tvalue.lensfile);
                lensfile      = lensfile.substr(0, lensfile.size() - 4);
                lensfile      = lensfile.substr(lensfile.find('.') + 1);
                lensfile      = lensfile.substr(0, lensfile.size() - 2);
                tvalue.approx_focallength = std::atof(lensfile.c_str());
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Film " + type);
    }
}

// Parse Texture
static inline void parse_texture(
    pbrt_stream& stream, const string& type, pbrt_texture& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "constant") {
        auto tvalue = pbrt_constant_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "value") {
                parse_param(stream, ptype, tvalue.value);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "bilerp") {
        auto tvalue = pbrt_bilerp_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "checkerboard") {
        auto tvalue = pbrt_checkerboard_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "dots") {
        auto tvalue = pbrt_dots_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "imagemap") {
        auto tvalue = pbrt_imagemap_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_mix_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "tex1") {
                parse_param(stream, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(stream, ptype, tvalue.tex2);
            } else if (pname == "amount") {
                parse_param(stream, ptype, tvalue.amount);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "scale") {
        auto tvalue = pbrt_scale_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "tex1") {
                parse_param(stream, ptype, tvalue.tex1);
            } else if (pname == "tex2") {
                parse_param(stream, ptype, tvalue.tex2);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fbm") {
        auto tvalue = pbrt_fbm_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "octaves") {
                parse_param(stream, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(stream, ptype, tvalue.roughness);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "wrinkled") {
        auto tvalue = pbrt_wrinkled_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "octaves") {
                parse_param(stream, ptype, tvalue.octaves);
            } else if (pname == "roughness") {
                parse_param(stream, ptype, tvalue.roughness);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "windy") {
        auto tvalue = pbrt_windy_texture{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "") {
                // TODO: missing params
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "marble") {
        auto tvalue = pbrt_marble_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uv") {
        auto tvalue = pbrt_uv_texture{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Texture " + type);
    }
}

variant<pbrt_plastic_material, pbrt_metal_material, pbrt_glass_material>
approximate_fourier_material(const string& filename) {
    if (get_filename(filename) == "paint.bsdf") {
        auto plastic = pbrt_plastic_material{};
        plastic.Kd   = {0.6f, 0.6f, 0.6f};
        // plastic.Ks = {0.4f, 0.4f, 0.4f};
        plastic.Ks         = {1.0f, 1.0f, 1.0f};
        plastic.uroughness = 0.2f;
        plastic.vroughness = 0.2f;
        return plastic;
    } else if (get_filename(filename) == "ceramic.bsdf") {
        auto plastic = pbrt_plastic_material{};
        plastic.Kd   = {0.6f, 0.6f, 0.6f};
        // plastic.Ks = {0.1f, 0.1f, 0.1f};
        plastic.Ks         = {1.0f, 1.0f, 1.0f};
        plastic.uroughness = 0.025f;
        plastic.vroughness = 0.025f;
        return plastic;
    } else if (get_filename(filename) == "leather.bsdf") {
        auto plastic = pbrt_plastic_material{};
        plastic.Kd   = {0.6f, 0.57f, 0.48f};
        // plastic.Ks = {0.1f, 0.1f, 0.1f};
        plastic.Ks         = {1.0f, 1.0f, 1.0f};
        plastic.uroughness = 0.3f;
        plastic.vroughness = 0.3f;
        return plastic;
    } else if (get_filename(filename) == "coated_copper.bsdf") {
        auto metal       = pbrt_metal_material{};
        auto etak        = pbrt_get_element_etak("Cu");
        auto [eta, k]    = etak;
        metal.eta        = {eta.x, eta.y, eta.z};
        metal.k          = {k.x, k.y, k.z};
        metal.uroughness = 0.01f;
        metal.vroughness = 0.01f;
        return metal;
    } else if (get_filename(filename) == "roughglass_alpha_0.2.bsdf") {
        auto glass       = pbrt_glass_material{};
        glass.uroughness = 0.2f;
        glass.vroughness = 0.2f;
        glass.Kr         = {1, 1, 1};
        glass.Kt         = {1, 1, 1};
        return glass;
    } else if (get_filename(filename) == "roughgold_alpha_0.2.bsdf") {
        auto metal       = pbrt_metal_material{};
        auto etak        = pbrt_get_element_etak("Au");
        auto [eta, k]    = etak;
        metal.eta        = {eta.x, eta.y, eta.z};
        metal.k          = {k.x, k.y, k.z};
        metal.uroughness = 0.2f;
        metal.vroughness = 0.2f;
        return metal;
    } else {
        throw io_error("unknown pbrt bsdf filename " + filename);
    }
}

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
static inline pair<vec3f, vec3f> parse_subsurface(const string& name) {
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
        {"Lowfat Milk",
            {{0.89187, 1.5136, 2.532}, {0.002875, 0.00575, 0.0115}}},
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
        {"Clorox", {{0.0024035, 0.0031373, 0.003991},
                       {0.0033542, 0.014892, 0.026297}}},
        {"Apple Juice", {{0.00013612, 0.00015836, 0.000227},
                            {0.012957, 0.023741, 0.052184}}},
        {"Cranberry Juice", {{0.00010402, 0.00011646, 7.8139e-05},
                                {0.039437, 0.094223, 0.12426}}},
        {"Grape Juice", {{5.382e-05, 0, 0}, {0.10404, 0.23958, 0.29325}}},
        {"Ruby Grapefruit Juice",
            {{0.011002, 0.010927, 0.011036}, {0.085867, 0.18314, 0.25262}}},
        {"White Grapefruit Juice",
            {{0.22826, 0.23998, 0.32748}, {0.0138, 0.018831, 0.056781}}},
        {"Shampoo", {{0.0007176, 0.0008303, 0.0009016},
                        {0.014107, 0.045693, 0.061717}}},
        {"Strawberry Shampoo", {{0.00015671, 0.00015947, 1.518e-05},
                                   {0.01449, 0.05796, 0.075823}}},
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
        {"Sugar Powder", {{0.00022272, 0.00025513, 0.000271},
                             {0.012638, 0.031051, 0.050124}}},
        {"Suisse Mocha Powder",
            {{2.7979, 3.5452, 4.3365}, {17.502, 27.004, 35.433}}},
        {"Pacific Ocean Surface Water", {{0.0001764, 0.00032095, 0.00019617},
                                            {0.031845, 0.031324, 0.030147}}},
    };
    return params.at(name);
}

// Get typename
static inline void parse_typeparam(pbrt_stream& stream, string& value) {
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
    if (value == "") throw io_error("type not found");
    restore_stream_position(stream);
}

// Parse param and resolve constant textures
static inline void parse_texture(pbrt_stream& stream, const string& ptype,
    pbrt_textured<pbrt_spectrum3f>&               value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
    parse_param(stream, ptype, value);
    if (value.texture == "") return;
    if (constant_values.find(value.texture) == constant_values.end()) return;
    value.value   = constant_values.at(value.texture);
    value.texture = "";
}
static inline void parse_texture(pbrt_stream& stream, const string& ptype,
    pbrt_textured<float>&                         value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
    parse_param(stream, ptype, value);
    if (value.texture == "") return;
    if (constant_values.find(value.texture) == constant_values.end()) return;
    auto col      = constant_values.at(value.texture);
    value.value   = (col.x + col.y + col.z) / 3;
    value.texture = "";
}

// Parse Material
static inline void parse_material(pbrt_stream& stream, const string& type,
    pbrt_material&                                value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
    auto pname = ""s, ptype = ""s;
    if (type == "matte") {
        auto tvalue = pbrt_matte_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "sigma") {
                parse_param(stream, ptype, tvalue.sigma);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mirror") {
        auto tvalue = pbrt_mirror_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kr") {
                parse_texture(stream, ptype, tvalue.Kr, constant_values);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "plastic") {
        auto tvalue = pbrt_plastic_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "Ks") {
                parse_texture(stream, ptype, tvalue.Ks, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "metal") {
        auto tvalue = pbrt_metal_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "k") {
                parse_texture(stream, ptype, tvalue.k, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "glass") {
        auto tvalue = pbrt_glass_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kr") {
                parse_texture(stream, ptype, tvalue.Kr, constant_values);
            } else if (pname == "Kt") {
                parse_texture(stream, ptype, tvalue.Kt, constant_values);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "translucent") {
        auto tvalue = pbrt_translucent_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "Ks") {
                parse_texture(stream, ptype, tvalue.Ks, constant_values);
            } else if (pname == "reflect") {
                parse_texture(stream, ptype, tvalue.reflect, constant_values);
            } else if (pname == "transmit") {
                parse_texture(stream, ptype, tvalue.transmit, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "uber") {
        auto tvalue = pbrt_uber_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "Ks") {
                parse_texture(stream, ptype, tvalue.Ks, constant_values);
            } else if (pname == "Kr") {
                parse_texture(stream, ptype, tvalue.Kr, constant_values);
            } else if (pname == "Kt") {
                parse_texture(stream, ptype, tvalue.Kt, constant_values);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "opacity") {
                parse_texture(stream, ptype, tvalue.opacity, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disney") {
        auto tvalue = pbrt_disney_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "color") {
                parse_texture(stream, ptype, tvalue.color, constant_values);
            } else if (pname == "anisotropic") {
                parse_texture(
                    stream, ptype, tvalue.anisotropic, constant_values);
            } else if (pname == "clearcoat") {
                parse_texture(stream, ptype, tvalue.clearcoat, constant_values);
            } else if (pname == "clearcoatgloss") {
                parse_texture(
                    stream, ptype, tvalue.clearcoatgloss, constant_values);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "metallic") {
                parse_texture(stream, ptype, tvalue.metallic, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "scatterdistance") {
                parse_texture(
                    stream, ptype, tvalue.scatterdistance, constant_values);
            } else if (pname == "sheen") {
                parse_texture(stream, ptype, tvalue.sheen, constant_values);
            } else if (pname == "sheentint") {
                parse_texture(stream, ptype, tvalue.sheentint, constant_values);
            } else if (pname == "spectrans") {
                parse_texture(stream, ptype, tvalue.spectrans, constant_values);
            } else if (pname == "thin") {
                parse_param(stream, ptype, tvalue.thin);
            } else if (pname == "difftrans") {
                parse_texture(stream, ptype, tvalue.difftrans, constant_values);
            } else if (pname == "flatness") {
                parse_texture(stream, ptype, tvalue.flatness, constant_values);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hair") {
        auto tvalue = pbrt_hair_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "color") {
                parse_texture(stream, ptype, tvalue.color, constant_values);
            } else if (pname == "sigma_a") {
                parse_texture(stream, ptype, tvalue.sigma_a, constant_values);
            } else if (pname == "eumelanin") {
                parse_texture(stream, ptype, tvalue.eumelanin, constant_values);
            } else if (pname == "pheomelanin") {
                parse_texture(
                    stream, ptype, tvalue.pheomelanin, constant_values);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "beta_m") {
                parse_texture(stream, ptype, tvalue.beta_m, constant_values);
            } else if (pname == "beta_n") {
                parse_texture(stream, ptype, tvalue.beta_n, constant_values);
            } else if (pname == "alpha") {
                parse_texture(stream, ptype, tvalue.alpha, constant_values);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "kdsubsurface") {
        auto tvalue = pbrt_kdsubsurface_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "Kr") {
                parse_texture(stream, ptype, tvalue.Kr, constant_values);
            } else if (pname == "Kt") {
                parse_texture(stream, ptype, tvalue.Kt, constant_values);
            } else if (pname == "mfp") {
                parse_texture(stream, ptype, tvalue.mfp, constant_values);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "mix") {
        auto tvalue = pbrt_mix_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "amount") {
                parse_texture(stream, ptype, tvalue.amount, constant_values);
            } else if (pname == "namedmaterial1") {
                parse_param(stream, ptype, tvalue.namedmaterial1);
            } else if (pname == "namedmaterial2") {
                parse_param(stream, ptype, tvalue.namedmaterial2);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "fourier") {
        auto tvalue = pbrt_fourier_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "bsdffile") {
                parse_param(stream, ptype, tvalue.bsdffile);
                tvalue.approx = approximate_fourier_material(tvalue.bsdffile);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "substrate") {
        auto tvalue = pbrt_substrate_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "Kd") {
                parse_texture(stream, ptype, tvalue.Kd, constant_values);
            } else if (pname == "Ks") {
                parse_texture(stream, ptype, tvalue.Ks, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "subsurface") {
        auto tvalue = pbrt_subsurface_material{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "name") {
                parse_param(stream, ptype, tvalue.name);
                auto params    = parse_subsurface(tvalue.name);
                tvalue.sigma_a = {
                    params.first.x, params.first.y, params.first.z};
                tvalue.sigma_prime_s = {
                    params.second.x, params.second.y, params.second.z};
            } else if (pname == "sigma_a") {
                parse_texture(stream, ptype, tvalue.sigma_a, constant_values);
            } else if (pname == "sigma_prime_s") {
                parse_texture(
                    stream, ptype, tvalue.sigma_prime_s, constant_values);
            } else if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "eta") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "index") {
                parse_texture(stream, ptype, tvalue.eta, constant_values);
            } else if (pname == "Kr") {
                parse_texture(stream, ptype, tvalue.Kr, constant_values);
            } else if (pname == "Kt") {
                parse_texture(stream, ptype, tvalue.Kt, constant_values);
            } else if (pname == "roughness") {
                pbrt_textured<float> roughness = 0.01f;
                parse_param(stream, ptype, roughness);
                tvalue.uroughness = roughness;
                tvalue.vroughness = roughness;
            } else if (pname == "uroughness") {
                parse_texture(
                    stream, ptype, tvalue.uroughness, constant_values);
            } else if (pname == "vroughness") {
                parse_texture(
                    stream, ptype, tvalue.vroughness, constant_values);
            } else if (pname == "remaproughness") {
                parse_param(stream, ptype, tvalue.remaproughness);
            } else if (pname == "bumpmap") {
                parse_texture(stream, ptype, tvalue.bumpmap, constant_values);
            } else if (pname == "type") {
                auto ttype = ""s;
                parse_param(stream, ptype, ttype);
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Material " + type);
    }
}

// Parse Shape
static inline void parse_shape(
    pbrt_stream& stream, const string& type, pbrt_shape& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "trianglemesh") {
        auto tvalue = pbrt_trianglemesh_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "plymesh") {
        auto tvalue = pbrt_plymesh_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "curve") {
        auto tvalue = pbrt_curve_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "loopsubdiv") {
        auto tvalue = pbrt_loopsubdiv_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "nurbs") {
        auto tvalue = pbrt_nurbs_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "sphere") {
        auto tvalue = pbrt_sphere_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "disk") {
        auto tvalue = pbrt_disk_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cone") {
        auto tvalue = pbrt_cone_shape{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "radius") {
                parse_param(stream, ptype, tvalue.radius);
            } else if (pname == "height") {
                parse_param(stream, ptype, tvalue.height);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "cylinder") {
        auto tvalue = pbrt_cylinder_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "hyperboloid") {
        auto tvalue = pbrt_hyperboloid_shape{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "p1") {
                parse_param(stream, ptype, tvalue.p1);
            } else if (pname == "p2") {
                parse_param(stream, ptype, tvalue.p2);
            } else if (pname == "phimax") {
                parse_param(stream, ptype, tvalue.phimax);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "paraboloid") {
        auto tvalue = pbrt_paraboloid_shape{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heightfield") {
        auto tvalue = pbrt_heightfield_shape{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "nu") {
                parse_param(stream, ptype, tvalue.nu);
            } else if (pname == "nv") {
                parse_param(stream, ptype, tvalue.nv);
            } else if (pname == "Pz") {
                parse_param(stream, ptype, tvalue.Pz);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Shape " + type);
    }
}

// Parse AreaLightSource
static inline void parse_arealight(
    pbrt_stream& stream, const string& type, pbrt_arealight& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "diffuse") {
        auto tvalue = pbrt_diffuse_arealight{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Film " + type);
    }
}

// Parse LightSource
static inline void parse_light(
    pbrt_stream& stream, const string& type, pbrt_light& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "distant") {
        auto tvalue = pbrt_distant_light{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "goniometric") {
        auto tvalue = pbrt_goniometric_light{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(stream, ptype, tvalue.I);
            } else if (pname == "mapname") {
                parse_param(stream, ptype, tvalue.mapname);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "infinite") {
        auto tvalue = pbrt_infinite_light{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "distant") {
        auto tvalue = pbrt_distant_light{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "L") {
                parse_param(stream, ptype, tvalue.L);
            } else if (pname == "from") {
                parse_param(stream, ptype, tvalue.from);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "projection") {
        auto tvalue = pbrt_projection_light{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "spot") {
        auto tvalue = pbrt_spot_light{};
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
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "point") {
        auto tvalue = pbrt_point_light{};
        while (is_param(stream)) {
            parse_nametype(stream, pname, ptype);
            if (pname == "scale") {
                parse_param(stream, ptype, tvalue.scale);
            } else if (pname == "I") {
                parse_param(stream, ptype, tvalue.I);
            } else if (pname == "from") {
                parse_param(stream, ptype, tvalue.from);
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown LightSource " + type);
    }
}

// Parse Medium
static inline void parse_pbrt_medium(
    pbrt_stream& stream, const string& type, pbrt_medium& value) {
    auto pname = ""s, ptype = ""s;
    if (type == "homogeneous") {
        auto tvalue = pbrt_homogeneous_medium{};
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
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else if (type == "heterogeneous") {
        auto tvalue = pbrt_heterogeneous_medium{};
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
                if (ttype != type) throw io_error("inconsistent types");
            } else {
                throw io_error("unknown parameter " + pname);
            }
        }
        value = tvalue;
    } else {
        throw io_error("unknown Medium " + type);
    }
}

// Load pbrt scene
template <typename Callbacks>
inline void load_pbrt(
    const string& filename, Callbacks& cb, const pbrt_params& params) {
    // start laoding files
    auto streams = vector<pbrt_stream>{};
    init_stream(streams);
    load_stream(filename, streams);

    // parsing stack
    auto stack    = vector<pbrt_context>{{}};
    auto object   = pbrt_object{};
    auto coordsys = unordered_map<string, pair<affine3f, affine3f>>{};

    // helpders
    auto set_transform = [](pbrt_context& ctx, const mat4f& xform) {
        if (ctx.active_transform_start) ctx.transform_start = (affine3f)xform;
        if (ctx.active_transform_end) ctx.transform_end = (affine3f)xform;
    };
    auto concat_transform = [](pbrt_context& ctx, const mat4f& xform) {
        if (ctx.active_transform_start) ctx.transform_start *= (affine3f)xform;
        if (ctx.active_transform_end) ctx.transform_end *= (affine3f)xform;
    };

    // constant values
    unordered_map<string, pbrt_spectrum3f> constant_values = {};

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
            if (stack.size() != 1) throw io_error("bad stack");
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
            stack.push_back(stack.back());
            cb.begin_object(object, stack.back());
        } else if (cmd == "ObjectEnd") {
            cb.end_object(object, stack.back());
            stack.pop_back();
            object = {};
        } else if (cmd == "ObjectInstance") {
            auto value = pbrt_object{};
            parse_value(stream, value.name);
            cb.object_instance(value, stack.back());
        } else if (cmd == "ActiveTransform") {
            auto value = ""s;
            parse_command(stream, value);
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
                throw io_error("bad active transform");
            }
        } else if (cmd == "Transform") {
            auto xf = identity_mat4f;
            parse_param(stream, xf);
            set_transform(stack.back(), xf);
        } else if (cmd == "ConcatTransform") {
            auto xf = identity_mat4f;
            parse_param(stream, xf);
            concat_transform(stack.back(), xf);
        } else if (cmd == "Scale") {
            auto v = zero3f;
            parse_param(stream, v);
            concat_transform(stack.back(), (mat4f)make_scaling_frame(v));
        } else if (cmd == "Translate") {
            auto v = zero3f;
            parse_param(stream, v);
            concat_transform(stack.back(), (mat4f)make_translation_frame(v));
        } else if (cmd == "Rotate") {
            auto v = zero4f;
            parse_param(stream, v);
            concat_transform(stack.back(),
                (mat4f)make_rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
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
            concat_transform(stack.back(), (mat4f)inverse(frame));
            stack.back().last_lookat_distance = length(from - to);
            // stack.back().focus = length(m.x - m.y);
        } else if (cmd == "ReverseOrientation") {
            stack.back().reverse = !stack.back().reverse;
        } else if (cmd == "CoordinateSystem") {
            auto name = ""s;
            parse_value(stream, name);
            coordsys[name] = {
                stack.back().transform_start, stack.back().transform_end};
        } else if (cmd == "CoordSysTransform") {
            auto name = ""s;
            parse_value(stream, name);
            if (coordsys.find(name) != coordsys.end()) {
                stack.back().transform_start = coordsys.at(name).first;
                stack.back().transform_end   = coordsys.at(name).second;
            }
        } else if (cmd == "Integrator") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_integrator{};
            parse_integrator(stream, type, value);
            cb.integrator(value, stack.back());
        } else if (cmd == "Sampler") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_sampler{};
            parse_sampler(stream, type, value);
            cb.sampler(value, stack.back());
        } else if (cmd == "PixelFilter") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_filter{};
            parse_filter(stream, type, value);
            cb.filter(value, stack.back());
        } else if (cmd == "Film") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_film{};
            parse_film(stream, type, value);
            cb.film(value, stack.back());
        } else if (cmd == "Accelerator") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_accelerator{};
            parse_accelerator(stream, type, value);
            cb.accelerator(value, stack.back());
        } else if (cmd == "Camera") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_camera{};
            parse_camera(stream, type, value);
            cb.camera(value, stack.back());
        } else if (cmd == "Texture") {
            auto name = ""s, comptype = ""s, type = ""s;
            parse_value(stream, name);
            parse_value(stream, comptype);
            parse_value(stream, type);
            auto value = pbrt_texture{};
            parse_texture(stream, type, value);
            if (type == "constant") {
                constant_values[name] =
                    get<pbrt_constant_texture>(value).value.value;
            }
            cb.texture(value, name, stack.back());
        } else if (cmd == "Material") {
            static auto material_id = 0;
            auto        type        = ""s;
            parse_value(stream, type);
            if (type == "") {
                stack.back().material = "";
            } else {
                auto value = pbrt_material{};
                auto name  = "unnamed_material_" + to_string(material_id++);
                parse_material(stream, type, value, constant_values);
                stack.back().material = name;
                cb.material(value, name, stack.back());
            }
        } else if (cmd == "MakeNamedMaterial") {
            auto name = ""s, type = ""s;
            parse_value(stream, name);
            parse_typeparam(stream, type);
            auto value = pbrt_material{};
            parse_material(stream, type, value, constant_values);
            cb.material(value, name, stack.back());
        } else if (cmd == "NamedMaterial") {
            auto name = ""s;
            parse_value(stream, name);
            stack.back().material = name;
        } else if (cmd == "Shape") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_shape{};
            parse_shape(stream, type, value);
            cb.shape(value, stack.back());
        } else if (cmd == "AreaLightSource") {
            auto type = ""s;
            parse_value(stream, type);
            static auto material_id = 0;
            auto        name  = "unnamed_arealight_" + to_string(material_id++);
            auto        value = pbrt_arealight{};
            parse_arealight(stream, type, value);
            stack.back().arealight = name;
            cb.arealight(value, name, stack.back());
        } else if (cmd == "LightSource") {
            auto type = ""s;
            parse_value(stream, type);
            auto value = pbrt_light{};
            parse_light(stream, type, value);
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
            load_stream(get_dirname(filename) + inputname, streams);
        } else {
            throw io_error("unknown command " + cmd);
        }
        // go to next file if needed
        skip_whitespace_or_comment_to_next_file(streams);
    }
}

}  // namespace yocto

#endif
