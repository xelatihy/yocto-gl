# Yocto/Trace: Path tracing

Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.

## Rendering features

Yocto/Trace uses the [Yocto/Scene](yocto_scene.md) model to represent
scene data. Consult the scene documentation for a detailed description
of the scene data model. Yocto/Trace supports rendering of shapes
represented as indexed meshes of points, lines, triangles and quads.
The renderer use instancing throughout to gain scalability with large
environment, while still remaining easy to use. Materials are modeled
similarly to the Disney's BSDF with integrated subsurface scattering.
Lights are not specified explicitly but derived by the renderer
from emissive surface and environment maps in the scene.

Yocto/Trace implements a fast but simple path tracer. The algorithm
uses relatively advanced sampling features to reduce the numbers of
rays used to render an image. At the same time, we focus on general
solutions that work well in all cases, instead of tuning the algorithm
to work well in some specific cases.

Yocto/Trace is a progressively path tracer that computes images
a few samples at a time. This is true for both offline and online
rendering. This design allows for partial images to be used to
give feedback before the render is finished. In the future, this
might allow for more rendering features like adaptive rendering
or stop-and-resume renders.

## Rendering a scene

To render a scene, first initialize the scene bvh and lights, with
`init_bvh(scene, params, progress)` and
`init_lights(scene, params, progress)` and then call
`trace_image(scene, camera, params, progress, image_progress)`.
In these functions, `params` are thee rendering options used to
render the scene, `progress` is a callback function used to
report rendering progress and `image_progress` is a callback
function used to report partial images during rendering.
Both callback functions are optional.

```cpp
auto scene = new scene_model{...};        // initialize scene
auto params = trace_params{};             // default params
auto progress = [](const string& message, // progress callback
      int current, int total) {
  print_info(message, current, total);
};
init_bvh(scene, params, progress);        // init bvh
init_lights(scene, params, progress);     // init lights
trace_image(scene, params, progress);     // render image
auto improgress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  if(sample % 16 == 0)                    // save image every 16 samples
    save_image(filename, render);
};
trace_image(scene, params, progress, improgress); // render image
```

## Rendering options

## Experimental async rendering

## Experimental denoising support

// Type of tracing algorithm
enum struct trace*sampler_type {
path, // path tracing
naive, // naive path tracing
eyelight, // eyelight rendering
falsecolor, // false color rendering
albedo, // renders the (approximate) albedo of objects for denoising
normal, // renders the normals of objects for denoising
};
// Type of false color visualization
enum struct trace_falsecolor_type {
// clang-format off
position, normal, frontfacing, gnormal, gfrontfacing, texcoord, color,
emission, diffuse, specular, coat, metal, transmission, translucency,
refraction, roughness, opacity, ior, instance, element, highlight
// clang-format on
};
// Strategy used to build the bvh
enum struct trace_bvh_type {
default*,
highquality,
middle,
balanced,
#ifdef YOCTO_EMBREE
embree_default,
embree_highquality,
embree_compact // only for copy interface
#endif
};

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Options for trace functions
struct trace*params {
int resolution = 1280;
trace_sampler_type sampler = trace_sampler_type::path;
trace_falsecolor_type falsecolor = trace_falsecolor_type::diffuse;
int samples = 512;
int bounces = 8;
float clamp = 100;
bool nocaustics = false;
bool envhidden = false;
bool tentfilter = false;
uint64_t seed = trace_default_seed;
trace_bvh_type bvh = trace_bvh_type::default*;
bool noparallel = false;
int pratio = 8;
float exposure = 0;
};

const auto trace_sampler_names = std::vector<std::string>{
"path", "naive", "eyelight", "falsecolor", "dalbedo", "dnormal"};

const auto trace_falsecolor_names = vector<string>{"position", "normal",
"frontfacing", "gnormal", "gfrontfacing", "texcoord", "color", "emission",
"diffuse", "specular", "coat", "metal", "transmission", "translucency",
"refraction", "roughness", "opacity", "ior", "instance", "element",
"highlight"};
const auto bvh_names = vector<string>{
"default", "highquality", "middle", "balanced",
#ifdef YOCTO_EMBREE
"embree-default", "embree-highquality", "embree-compact"
#endif
};

// Progress report callback
using progress_callback =
function<void(const string& message, int current, int total)>;
// Callback used to report partially computed image
using image_callback =
function<void(const image<vec4f>& render, int current, int total)>;

// Initialize lights.
void init_lights(scene_model\* scene, progress_callback progress_cb = {});

// Build the bvh acceleration structure.
void init_bvh(scene_model\* scene, const trace_params& params,
progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(scene_model* scene,
const vector<scene_instance*>& updated_objects,
const vector<scene_shape\*>& updated_shapes, const trace_params& params);

// Progressively computes an image.
image<vec4f> trace_image(const scene_model* scene, const scene_camera* camera,
const trace_params& params, progress_callback progress_cb = {},
image_callback image_cb = {});

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// [experimental] Asynchronous state
struct trace_state {
image<vec4f> render = {};
image<vec4f> accumulation = {};
image<int> samples = {};
image<rng_state> rngs = {};
future<void> worker = {}; // async
atomic<bool> stop = {}; // async
};

// [experimental] Callback used to report partially computed image
using async_callback = function<void(
const image<vec4f>& render, int current, int total, const vec2i& ij)>;

// [experimental] Asynchronous interface
struct trace_state;
void trace_start(trace_state* state, const scene_model* scene,
const scene_camera* camera, const trace_params& params,
progress_callback progress_cb = {}, image_callback image_cb = {},
async_callback async_cb = {});
void trace_stop(trace_state* state);
