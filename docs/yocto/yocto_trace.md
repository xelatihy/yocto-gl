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
  if(sample % 16 == 0)                    // every 16 samples
    save_image(filename, render);         // save image
};
trace_image(scene, params, progress, improgress); // render image
```

## Rendering options

The rendering process is configured with the `trace_params` settings.
`sampler` determines the algorithm used for rendering: `path` is the default
algorithm and probably the one you want to use, `naive` is a simpler path
tracing that may be used for testing, `eyelight` produces quick previews
of the screen geometry, `falsecolor` is a debug feature to view scenes
according to the `falsecolor` setting, and `albedo` and `normal` are
used to produce denoising buffers.

THe image resolution is set by `resolution` and measures the resolution
of the longest axis. `samples` is the number of per-pixel samples
used while rendering and is the only parameter used to control the
tradeoff between noise and speed. `bounces` is the maximum number of bounces
and should be high for scenes with glass and volumes, but otherwise a low
number would suffice.

The remaining parameters are approximation used to reduce noise, at the
expenses of bias. `clamp` remove high-energy fireflies. `nocaustics` removes
certain path that cause caustics. `tentfilter` apply a linear filter to the
image pixels. `envhidden` removes the environment map from the camera rays.

Finally, the `bvh` parameter controls the heuristic used to build the Bvh
and whether the Bvh uses Embree. Please see the description in Yocto/Scene.

`trace_sampler_names`, `trace_falsecolor_names` and `trace_bvh_names`
define string names for various enum values that can used for UIs or CLIs.

```cpp
// high quality rendering
auto hq_params = trace_params{};             // default params
hq_params.sampler    = trace_sampler_type::path; // path tracing
hq_params.resolution = 1280;                 // high-res render
hq_params.samples    = 1024;                 // high-sample count
hq_params.bounce     = 64;                   // high max bounces
// geometry previewing
auto pp_params = trace_params{};             // default params
pp_params.sampler    = trace_sampler_type::eyelight; // geometry preview
pp_params.resolution = 720;                  // medium-res render
pp_params.samples    = 16;                   // low-sample count
// scene debug viewing
auto db_params = trace_params{};             // default params
db_params.sampler    = trace_sampler_type::falsecolor; // debug false colors
db_params.falsecolor = trace_falsecolor_type::normal; // normal viewing
db_params.resolution = 720;                  // medium-res render
db_params.samples    = 16;                   // low-sample count
```

## Experimental async rendering

The render can run in asynchronous mode where the rendering process is
lunched and runs until completion. This mode is useful when for interactive
viewing or in modeling-while-rendering applications. The API is very minimal
and only controls the rendering process. Use `trace_start(...)` to start
the async renderer and `trace_stop(...)` to stop it. The renderer starts
by rendering a low resolution preview and then proceeds progressively.

The async renderer takes a `trace_state` struct that tracks the rendering
process and contains all data needed by the async renderer. Rendering progress
ifs given by three callbacks. The first two are the rendering callbacks
defined for offline rendeirng, that return progress report and an image buffer
for each sample. In async mode, a further callback is called after each
pixel is rendered.

During rendering, no scenes changes are allowed, and changes to bvh and lights
are not tracked. This is on purpose since it allows for a simple API while
retaining maximum speed. To changes the scene, first stop the render, than
apply scene changes, than update lights and bvh if needed, and finally restart
the renderer.
To update the BVH, use `update_bvh(scene, instances, shapes, params)` where
`instances` and `shapes` are the list of modified ids.

```cpp
auto scene = new scene_model{...};        // initialize scene
auto params = trace_params{};             // default params
init_bvh(scene, params, progress);        // init bvh
init_lights(scene, params, progress);     // init lights
auto imgprogress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  display_image(render);                  // display image
};
auto pxlprogress = [](const image<vec4f>& render, // pixel progress
      int sample, int samples, const vec2i& ij) {
  display_pixel(render, ij);              // display pixel
};
auto state = new trace_state{};           // allocate state
trace_start(state, scene, params, {},     // start async renderer
  imgprogress, pxlprogress);              // function returns immediately
run_app(...);                             // application runs normally here
trace_stop(state);                        // stop async renderer
modify_scene(...);                        // make scene changes
trace_start(state, scene, params, {},     // re-start async renderer
  imgprogress, pxlprogress);              // function returns immediately
```
