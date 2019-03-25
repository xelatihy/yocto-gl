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

#include <unordered_map>
#include <variant>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::function;
using std::unordered_map;
using std::variant;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// pbrt cameras
struct pbrt_camera_perspective {
    float fov              = 90;
    float frameaspectratio = -1;  // or computed from film
    float lensradius       = 0;
    float focaldistance    = 1e30;
    vec4f screenwindow     = {-1, -1, 1, 1};
    float shutteropen      = 0;
    float shutterclose     = 1;
};
struct pbrt_camera_orthographic {
    float frameaspectratio = -1;  // or computed from film
    float lensradius       = 0;
    float focaldistance    = 1e30;
    vec4f screenwindow     = {-1, -1, 1, 1};
    float shutteropen      = 0;
    float shutterclose     = 1;
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
struct pbrt_film {
    int    xresolution        = 640;
    int    yresolution        = 480;
    vec4f  cropwindow         = {0, 1, 0, 1};
    float  scale              = 1;
    float  maxsampleluminance = type_max<float>();
    float  diagonal           = 35;
    string filename           = "pbrt.exr";
};

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
    int                   maxdepth            = 5;
    vec4i                 pixelbounds         = {-1, -1, -1, -1};
    float                 rrthreshold         = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
};
struct pbrt_integrator_bdpt {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int                   maxdepth            = 5;
    vec4i                 integer             = {-1, -1, -1, -1};
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::power;
    bool                  visualizestrategies = false;
    bool                  visualizeweights    = false;
};
struct pbrt_integrator_directlighting {
    enum struct strategy_t { all, one };
    strategy_t strategy    = strategy_t::all;
    int        maxdepth    = 5;
    vec4i      pixelbounds = {-1, -1, -1, -1};
};
struct pbrt_integrator_mlt {
    int   maxdepth             = 5;
    int   bootstrapsamples     = 100000;
    int   chains               = 1000;
    int   mutationsperpixel    = 100;
    float largestepprobability = 0.3;
    float sigma                = 0.01;
};
struct pbrt_integrator_sppm {
    int   maxdepth            = 5;
    int   iterations          = 64;
    int   photonsperiteration = -1;
    int   imagewritefrequency = pow2(31);
    float radius              = 5;
};
struct pbrt_integrator_whitted {
    // TODO: missing from documentation
};
using pbrt_integrator = variant<pbrt_integrator_path, pbrt_integrator_bdpt,
    pbrt_integrator_directlighting, pbrt_integrator_mlt, pbrt_integrator_sppm,
    pbrt_integrator_whitted>;

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

// pbrt textures
struct pbrt_texture_constant {
    vec3f value = {1, 1, 1};
};
struct pbrt_texture_bilerp {
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
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_fbm {};
struct pbrt_texture_imagemap {
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
};
struct pbrt_texture_marble {};
struct pbrt_texture_mix {};
struct pbrt_texture_scale {};
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
struct pbrt_texture_windy {};
struct pbrt_texture_wrinkled {};
using pbrt_texture = variant<pbrt_texture_constant, pbrt_texture_bilerp,
    pbrt_texture_checkerboard, pbrt_texture_dots, pbrt_texture_fbm,
    pbrt_texture_imagemap, pbrt_texture_marble, pbrt_texture_mix,
    pbrt_texture_scale, pbrt_texture_uv, pbrt_texture_windy,
    pbrt_texture_wrinkled>;

// pbrt materials
struct pbrt_material_none {};
struct pbrt_material_matte {
    string       name      = "";
    vec3f        Kd        = {0.5f, 0.5f, 0.5f};
    pbrt_texture Kd_map    = {};
    float        sigma     = 0;
    pbrt_texture sigma_map = {};
};
struct pbrt_material_mirror {};
struct pbrt_material_plastic {};
struct pbrt_material_metal {};
struct pbrt_material_glass {};
struct pbrt_material_translucent {};
struct pbrt_material_uber {};
struct pbrt_material_disney {};
struct pbrt_material_fourier {};
struct pbrt_material_hair {};
struct pbrt_material_kdsubsurface {};
struct pbrt_material_mix {};
struct pbrt_material_substrate {};
struct pbrt_material_subsurface {};
using pbrt_material =
    variant<pbrt_material_none, pbrt_material_matte, pbrt_material_mirror,
        pbrt_material_plastic, pbrt_material_metal, pbrt_material_glass,
        pbrt_material_translucent, pbrt_material_uber, pbrt_material_disney,
        pbrt_material_fourier, pbrt_material_hair, pbrt_material_kdsubsurface,
        pbrt_material_mix, pbrt_material_substrate, pbrt_material_subsurface>;

// pbrt shapes
struct pbrt_shape_trianglemesh {
    vector<vec3i> indices = {};
    vector<vec3f> P       = {};
    vector<vec3f> N       = {};
    vector<vec3f> S       = {};
    vector<vec2f> uv      = {};
    // texture alpha
    // texture shadowalpha
};
struct pbrt_shape_plymesh {
    string filename = {};
    // texture alpha
    // texture shadowalpha
};
struct pbrt_shape_curve {
    enum struct type_t { flat, ribbon, cylinder };
    enum struct basis_t { bezier, bspline };
    vector<vec3f> P          = {};
    basis_t       basis      = basis_t::bezier;
    int           degree     = 3;
    type_t        type       = type_t::flat;
    vector<vec3f> N          = {};
    float         width      = 1;
    float         width0     = 1;
    float         width1     = 1;
    int           splitdepth = 3;
};
struct pbrt_shape_loopsubdiv {
    int           levels  = 3;
    vector<int>   indices = {};
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
    vec3f scale = {1, 1, 1};
    vec3f L     = {1, 1, 1};
    vec3f from{0, 0, 0};
    vec3f to = {0, 0, 1};
};
struct pbrt_light_goniometric {
    vec3f  scale   = {1, 1, 1};
    vec3f  I       = {1, 1, 1};
    string mapname = "";
};
struct pbrt_light_infinite {
    vec3f  scale   = {1, 1, 1};
    vec3f  L       = {1, 1, 1};
    int    samples = 1;
    string mapname = "";
};
struct pbrt_light_point {
    vec3f scale = {1, 1, 1};
    vec3f I     = {1, 1, 1};
    vec3f from{0, 0, 0};
};
struct pbrt_light_projection {
    vec3f  scale   = {1, 1, 1};
    vec3f  I       = {1, 1, 1};
    float  fov     = 45;
    string mapname = "";
};
struct pbrt_light_spot {
    vec3f scale          = {1, 1, 1};
    vec3f I              = {1, 1, 1};
    vec3f from           = {0, 0, 0};
    vec3f to             = {0, 0, 1};
    float coneangle      = 30;
    float conedeltaangle = 5;
};
using pbrt_light =
    variant<pbrt_light_distant, pbrt_light_goniometric, pbrt_light_infinite,
        pbrt_light_point, pbrt_light_projection, pbrt_light_spot>;

// pbrt area lights
struct pbrt_arealight_none {};
struct pbrt_arealight_diffuse {
    vec3f L        = {1, 1, 1};
    bool  twosided = false;
    int   samples  = 1;
};
using pbrt_arealight = variant<pbrt_arealight_none, pbrt_arealight_diffuse>;

// pbrt callbacks
struct pbrt_callbacks {};

// Load obj options
struct load_pbrt_options {
    bool geometry_only = false;
    bool flip_texcoord = true;
};

// Load obj scene
void load_pbrt(const string& filename, const pbrt_callbacks& cb,
    const load_pbrt_options& options = {});

// objio error
struct pbrtio_error : runtime_error {
    explicit pbrtio_error(const char* msg) : runtime_error{msg} {}
    explicit pbrtio_error(const std::string& msg) : runtime_error{msg} {}
};

}  // namespace yocto

#endif
