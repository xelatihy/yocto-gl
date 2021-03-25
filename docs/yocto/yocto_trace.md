# Yocto/Trace: Path tracing

Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.

## Rendering features

Yocto/Trace uses the [Yocto/Scene](yocto_scene.md) model to represent
scene data. Consult the scene documentation for a detailed description
of the scene data model.

Yocto/Trace implements a fast but simple path tracer. The algorithm
uses relatively advanced sampling features to reduce the numbers of
rays used to render an image. At the same time, we focus on general
solutions that work well in all cases, instead of tuning the algorithm
to work well in some specific cases.

Yocto/Trace is a progressive path tracer that computes images
a few samples at a time. This is true for both offline and online
rendering. This design allows for partial images to be used to
give feedback before the render is finished.

## Rendering a scene

To render a scene, first tesselate shapes for subdivs and displacement,
with `tesselate_shapes(scene, params)`, then call
`trace_image(scene, camera, params)`, where `params` are the rendering options.

```cpp
auto scene = scene_model{...};            // initialize scene
auto params = trace_params{};             // default params
tesselate_shapes(scene, params);          // tesselate shapes if needed
trace_image(scene, params);               // render image
```

## Rendering options

The rendering process is configured with the `trace_params` settings.
`sampler` determines the algorithm used for rendering: `path` is the default
algorithm and probably the one you want to use, `naive` is a simpler path
tracing that may be used for testing, `eyelight` produces quick previews
of the screen geometry, `falsecolor` is a debug feature to view scenes
according to the `falsecolor` setting.

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

Finally, `highqualitybvh` congtrols the BVH quality and `embreebvh` controls
whether to use Intel's Embree. Please see the description in
[Yocto/Bvh](yocto_bvh.md).

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

## Low-Level Rendering API

Yocto/Trace supports a low-level rendering API that is more flexible for
applications that need it. In this modality, you have to manage all state
manually.

Yocto/Trace internally uses a `bvh_scene` BVH for ray-scene intersection,
a `trace_lights` objects to store lighting information, and a `trace_state`
object to tracks the rendering process and contains all data needed to perform
progressive computation. Each object need to be initialized separately.

To render a scene, first tesselate shapes for subdivs and displacement,
with `tesselate_shapes(scene, params, progress)`, then initialize the scene
bvh and lights, with `make_bvh(scene, params)` and
`make_lights(scene, params)`, then the rendering state
with `make_state(state, scene)`.

Then, for each sample, call `trace_samples(state, scene, lights, params)`
and retrieve the computed image with `get_render(state)` or
`get_render(image, state)`. This interface can be useful to provide user
feedback by either saving or displaying partial images.

```cpp
auto scene = scene_model{...};              // initialize scene
auto params = trace_params{};               // default params
tesselate_shapes(scene, params);            // tesselate shapes if needed
auto bvh = make_bvh(scene, params);         // init bvh
auto lights = make_lights(scene, params);   // init lights
auto state = make_state(scene, params);     // init state
for(auto sample : range(params.samples)) {  // for each sample
  trace_samples(state, scene, camera, bvh,  // render sample
                lights, params);
  process_image(get_render(state));          // get image computed so far
};
```

## Denoising with Intel's Open Image Denoise

We support denoising of rendered images in the low-level interface.
Just call `get_denoised(...)` instead of `get_render(...)` to get a denoised image.
Alternatively, you can call `get_albedo(state)` or `get_normal(state)` to get
denoising buffers and either run the denoiser in a different process, or
call `denoise_render(render, albedo, normal)` to denoise the image.
To denoise within Yocto/GL, the library should be compiled with OIDN support by
setting the `YOCTO_DENOISE` compile flag and linking to OIDN's libraries.

```cpp
auto scene = scene_model{...};              // initialize scene
auto params = trace_params{};               // default params
tesselate_shapes(scene, params);            // tesselate shapes if needed
auto bvh = make_bvh(scene, params);         // init bvh
auto lights = make_lights(scene, params);   // init lights
auto state = make_state(scene, params);     // init state
for(auto sample : range(params.samples)) {  // for each sample
  trace_samples(state, scene, camera, bvh,  // render sample
                lights, params);
};
auto denoised = get_denoised(state);        // get denoised final image
// alternative interface
auto render = get_render(state);            // get final image
auto albedo = get_albedo(state);            // get denoising buffers
auto normal = get_normal(state);
// run denoiser here or save buffers and run elsewhere
auto denoised2 = denoise_render(render, albedo, normal);
```
