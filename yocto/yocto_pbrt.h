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

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::function;
using std::unordered_map;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// pbrt cameras
struct pbrt_camera {
    struct perspective_t {
        float fov              = 90;
        float frameaspectratio = -1;  // or computed from film
        float lensradius       = 0;
        float focaldistance    = 1e30;
        vec4f screenwindow     = {-1, -1, 1, 1};
        float shutteropen      = 0;
        float shutterclose     = 1;
    };
    struct orthographic_t {
        float frameaspectratio = -1;  // or computed from film
        float lensradius       = 0;
        float focaldistance    = 1e30;
        vec4f screenwindow     = {-1, -1, 1, 1};
        float shutteropen      = 0;
        float shutterclose     = 1;
    };
    struct environment_t {
        float shutteropen  = 0;
        float shutterclose = 1;
    };
    struct realistic_t {
        string lensfile         = "";
        float  aperturediameter = 1;
        float  focusdistance    = 10;
        bool   simpleweighting  = true;
        float  shutteropen      = 0;
        float  shutterclose     = 1;
    };
    enum struct type_t { perspective, orthographic, environment, realistic };
    type_t         type         = type_t::perspective;
    perspective_t  perspective  = {};
    orthographic_t orthographic = {};
    environment_t  environment  = {};
    realistic_t    realistic    = {};
};

// pbrt samplers
struct pbrt_sampler {
    struct random_t {
        int pixelsamples = 16;
    };
    struct halton_t {
        int pixelsamples = 16;
    };
    struct sobol_t {
        int pixelsamples = 16;
    };
    struct zerotwosequence_t {
        int pixelsamples = 16;
    };
    struct maxmindist_t {
        int pixelsamples = 16;
    };
    struct stratified_t {
        bool jitter   = true;
        int  xsamples = 2;
        int  ysamples = 2;
    };
    enum struct type_t {
        random,
        halton,
        sobol,
        zerotwosequence,
        maxmindist,
        stratified
    };
    type_t            type            = type_t::random;
    random_t          random          = {};
    halton_t          halton          = {};
    sobol_t           sobol           = {};
    zerotwosequence_t zerotwosequence = {};
    maxmindist_t      maxmindist      = {};
    stratified_t      stratified      = {};
};

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
struct pbrt_filter {
    struct box_t {
        float xwidth = 0.5f;
        float ywidth = 0.5f;
    };
    struct gaussian_t {
        float xwidth = 2;
        float ywidth = 2;
        float alpha  = 2;
    };
    struct mitchell_t {
        float xwidth = 2;
        float ywidth = 2;
        float B      = 1.0f / 3.0f;
        float C      = 1.0f / 3.0f;
    };
    struct sinc_t {
        float xwidth = 4;
        float ywidth = 4;
        float tau    = 3;
    };
    struct triangle_t {
        float xwidth = 2;
        float ywidth = 2;
    };
    enum struct type_t { none, box, guassian, mitchell, sinc, triangle };
    type_t     type     = type_t::none;
    box_t      box      = {};
    gaussian_t guassian = {};
    mitchell_t mitchell = {};
    sinc_t     sinc     = {};
    triangle_t triangle = {};
};

// pbrt integrators
struct pbrt_integrator {
    struct path_t {
        enum struct lightsamplestrategy_t { uniform, power, spatial };
        int                   maxdepth    = 5;
        vec4i                 pixelbounds = {-1, -1, -1, -1};
        float                 rrthreshold = 1;
        lightsamplestrategy_t lightsamplestrategy =
            lightsamplestrategy_t::spatial;
    };
    struct bdpt_t {
        enum struct lightsamplestrategy_t { uniform, power, spatial };
        int                   maxdepth = 5;
        vec4i                 integer  = {-1, -1, -1, -1};
        lightsamplestrategy_t lightsamplestrategy =
            lightsamplestrategy_t::power;
        bool visualizestrategies = false;
        bool visualizeweights    = false;
    };
    struct directlighting_t {
        enum struct strategy_t { all, one };
        strategy_t strategy    = strategy_t::all;
        int        maxdepth    = 5;
        vec4i      pixelbounds = {-1, -1 -, 1 -, 1};
    };
    struct mlt_t {
        int   maxdepth             = 5;
        int   bootstrapsamples     = 100000;
        int   chains               = 1000;
        int   mutationsperpixel    = 100;
        float largestepprobability = 0.3;
        float sigma                = 0.01;
    };
    struct sppm_t {
        int   maxdepth            = 5;
        int   iterations          = 64;
        int   photonsperiteration = -1;
        int   imagewritefrequency = pow2(31);
        float radius              = 5;
    };
    struct whitted_t {
        // TODO: missing from documentation
    };
    enum struct type_t { path, bdpt, directlighting, mlt, sppm, whitted };
    type_t           type           = type_t::path;
    path_t           path           = {};
    bdpt_t           bdpt           = {};
    directlighting_t directlighting = {};
    mlt_t            mlt            = {};
    sppm_t           sppm           = {};
    whitted_t        whitted        = {};
};

// pbrt accellerators
struct pbrt_accelerator {
    struct bvh_t {
        enum struct splitmethod_t { sah, equal, middle, hlbvh };
        int           maxnodeprims = 4;
        splitmethod_t splitmethod  = splitmethod_t::sah;
    };
    struct kdtree_t {
        int   intersectcost = 80;
        int   traversalcost = 1;
        float emptybonus    = 0.2;
        int   maxprims      = 1;
        int   maxdepth      = -1;
    };
    enum struct type_t { bvh, kdtree };
    type_t   type   = type_t::bvh;
    bvh_t    bvh    = {};
    kdtree_t kdtree = {};
};

// pbrt textures
struct pbrt_texture {
    
};

// pbrt materials
struct pbrt_material {
struct none_t {
};
struct matte_t {
};
struct mirror_t {
};
struct plastic_t {
};
struct metal_t {
};
struct glass_t {
};
struct translucent_t {
};
struct uber_t {
};
struct disney_t {
};
struct fourier_t {
};
struct hair_t {
};
struct kdsubsurface_t {
};
struct mix_t {
};
struct substrate_t {
};
struct subsurface_t {
};
    enum type_t {
none,
matte,
mirror,
plastic,
metal,
glass,
translucent,
uber,
disney,
fourier,
hair,
kdsubsurface,
mix,
substrate,
subsurface,
    };

    string name = "";
    type_t type = type_t::matte;
none_t none = {};
matte_t matte = {};
mirror_t mirror = {};
plastic_t plastic = {};
metal_t metal = {};
glass_t glass = {};
translucent_t translucent = {};
uber_t uber = {};
disney_t disney = {};
fourier_t fourier = {};
hair_t hair = {};
kdsubsurface_t kdsubsurface = {};
mix_t mix = {};
substrate_t substrate = {};
subsurface_t subsurface = {};
};

// pbrt shapes
struct pbrt_shape {
    struct trianglemesh_t {
        vector<vec3i> indices = {};
        vector<vec3f> P       = {};
        vector<vec3f> N       = {};
        vector<vec3f> S       = {};
        vector<vec2f> uv      = {};
        // texture alpha
        // texture shadowalpha
    };
    struct plymesh_t {
        string filename = {};
        // texture alpha
        // texture shadowalpha
    };
    struct curve_t {
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
    struct loopsubdiv_t {
        int           levels  = 3;
        vector<int>   indices = {};
        vector<vec3f> P       = {};
    };
    struct nurbs_t {
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
    struct sphere_t {
        float radius = 1;
        float zmin   = -radius;
        float zmax   = radius;
        float phimax = 360;
    };
    struct disk_t {
        float height      = 0;
        float radius      = 1;
        float innerradius = 0;
        float phimax      = 360;
    };
    struct cone_t {
        float radius = 1;
        float height = 1;
        float phimax = 360;
    };
    struct cylinder_t {
        float radius = 1;
        float zmin   = -1;
        float zmax   = 1;
        float phimax = 360;
    };
    struct hyperboloid_t {
        vec3f p1     = {0, 0, 0};
        vec3f p2     = {1, 1, 1};
        float phimax = 360;
    };
    struct paraboloid_t {
        float radius = 1;
        float zmin   = 0;
        float zmax   = 1;
        float phimax = 360;
    };
    struct heightfield_t {
        int           nu = 0;
        int           nv = 0;
        vector<float> Pz = {};
    };
    enum struct type_t {
        trianglemesh,
        plymesh,
        curve,
        loopsubdiv,
        nurbs,
        sphere,
        disk,
        cone,
        cylinder,
        hyperboloid,
        paraboloid,
        heightfield
    };
    type_t         type         = type_t::trianglemesh;
    trianglemesh_t trianglemesh = {};
    plymesh_t      plymesh      = {};
    curve_t        curve        = {};
    loopsubdiv_t   loopsubdiv   = {};
    nurbs_t        nurbs        = {};
    sphere_t       sphere       = {};
    disk_t         disk         = {};
    cone_t         cone         = {};
    cylinder_t     cylinder     = {};
    hyperboloid_t  hyperboloid  = {};
    paraboloid_t   paraboloid   = {};
    heightfield_t  heightfield  = {};
};

// pbrt lights
struct pbrt_light {
    struct distant_t {
        vec3f scale = {1, 1, 1};
        vec3f L     = {1, 1, 1};
        vec3f from{0, 0, 0};
        vec3f to = {0, 0, 1};
    };
    struct goniometric_t {
        vec3f  scale   = {1, 1, 1};
        vec3f  I       = {1, 1, 1};
        string mapname = "";
    };
    struct infinite_t {
        vec3f  scale   = {1, 1, 1};
        vec3f  L       = {1, 1, 1};
        int    samples = 1;
        string mapname = "";
    };
    struct point_t {
        vec3f scale = {1, 1, 1};
        vec3f I     = {1, 1, 1};
        vec3f from{0, 0, 0};
    };
    struct projection_t {
        vec3f  scale   = {1, 1, 1};
        vec3f  I       = {1, 1, 1};
        float  fov     = 45;
        string mapname = "";
    };
    struct spot_t {
        vec3f scale          = {1, 1, 1};
        vec3f I              = {1, 1, 1};
        vec3f from           = {0, 0, 0};
        vec3f to             = {0, 0, 1};
        float coneangle      = 30;
        float conedeltaangle = 5;
    };
    enum struct type_t {
        distant,
        goniometric,
        infinite,
        point,
        projection,
        spot
    };
    type_t        type        = type_t::distant;
    distant_t     distant     = {};
    goniometric_t goniometric = {};
    infinite_t    infinite    = {};
    point_t       point       = {};
    projection_t  projection  = {};
    spot_t        spot        = {};
};

// pbrt area lights
struct pbrt_arealight {
    struct diffuse_t {
        vec3f	L	= {1, 1, 1};
        bool	twosided =	false;
        int	samples	= 1;
    };
    enum struct type_t { none, diffuse };
    type_t type = type_t::none;
    diffuse_t diffuse = {};
};

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
