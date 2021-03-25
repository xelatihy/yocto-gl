# Yocto/Image: Image utilities

Yocto/Image is a collection of image utilities useful when writing rendering
algorithms. These include a simple image data structure, color conversion
utilities and tone mapping, and image resizing.
Yocto/Image is implemented in `yocto_image.h` and `yocto_image.cpp`, and
depends on `stb_image_resize.h`.

## Image representation

Images are represented as a simple struct called `color_image`, that stores
the image width and height, a boolean flat that indicates whether the colors
are in linear ot sRGB color space, and an array of pixels of `vec4f` type.

While images can be constructed by direct manipulation of their values,
th preferred way to construct images is by using the
`make_image(widgh, height, linear)` function, where colors are default
initialized to black.

Image pixels are stored in row-major order and be accessed with index
`j * width + i`. They can be conveniently accessed using `get_pixel(image,i,j)`
and `set_pixel(image,i,j,color)` function.

```cpp
auto image = make_image(512,512,false); // creates a 512x512 sRGB image
for(auto& pixel : image.pixels)         // iterate over pixels
  print_info(pixel);
auto color = get_pixel(image,128,128);  // pixel access
for(auto j=0; j<image.height; j++)      // iterate over image rows
  for(auto i=0; i<image.width; i++)     // iterate over row pixels
    print_info(get_pixel(image,i,j));
```

Images are sampled with `eval_image(img,uv,...)`, where `uv` are in
the [0,1]x[0,1] range. By default, image values are returned in the color
space of the image, but can be forced to to be linear.
By default, images are sampled with bilinear interpolation and tiling,
with nearest-neighbor interpolation and edge-clamp behavior available
as overrides.

```cpp
auto hdr = make_image(w,h,true);        // creates a linear image
auto ch0 = eval_image(hdr, {0.5,0.5});  // samples the image center
auto ch1 = eval_image(hdr, {0.5,0.5}, false, false);  // no interpolation
auto ch2 = eval_image(hdr, {0.5,0.5}, false, false, false);  // clamp to edge
auto ldr = make_image(w,h,false);       // creates a sRGB image
auto ch0 = eval_image(ldr, {0.5,0.5});  // samples in sRGB, returns as sRGB
auto ch1 = eval_image(ldr, {0.5,0.5}, true);  // treats sRGB values as linear
```

Image loading and saving is defined in [Yocto/SceneIO](yocto_sceneio.md).

## Image utilities

Images can be converted between linear RGB and sRGB color space using
`convert_image(image,linear)`.

```cpp
auto linear = make_image(w,h,true);   // initialize a linear image
auto srgb = convert_image(w,h,false); // convert to sRGB
```

HDR images can be tone mapped using `tonemap_image(hdr,exposure,filmic)`
that applies an exposure correction followed by an optional filmic tone curve.
Use `tonemap_image_mt(...)` to quickly tone map large images using
multiple threads.

```cpp
auto hdr = make_image(w,h,true);         // initialize am HDR image
auto ldr = tonemap_image(hdr, 2);        // tone mapping with exposure 2
auto flm = tonemap_image(hdr, 2, true);  // filmic tone mapping
```

Images can be color graded by applying a set of minimal color grading tools
using `colorgrade_image(image,params)`, in manner similar to
[Hable 2017](http://filmicworlds.com/blog/minimal-color-grading-tools/).
Color grading corrections can be applied on images that are either linear HDR
or sRGB encoded. The results is always an LDR image encoded in sRGB.
Use `colorgrade_image_mt(...)` to quickly tone map large images using
multiple threads.

Several color corrections are bundled together with their parameters
packed in a convenience structure `colorgrade_params`.
Color grading operations are applied in a fixed sequence and consist of the
following operations: exposure compensation, color tint, contrast in the log domain,
filmic curve, conversion to sRGB, S-shaped contrast, saturation,
and shadow/midtone/highlight correction.
Color tinting can be used to apply white balance by using
`compute_white_balance(img)` to determine the correct color.

```cpp
auto hdr = make_image(w,h,true);         // initialize an HDR image
auto params = colorgrade_params{};       // default color grading params
params.exposure = 1;                     // set desired params
params.logcontrast = 0.75;               // set desired params
params.tint = compute_white_balance(hdr);// apply white balance
auto ldr = colorgrade_image(hdr, params);// color grading
```

Images are can resized with `resize_image(image,w,h)`. Just like all other
functions, images are not resized in placed, but a new image is created.
Resizing works for both linear and 8bit images.

```cpp
auto img = make_image(...);              // initialize an image
auto res = resize_image(img, 512, 512);  // resizing to fixed size
auto asp = resize_image(img, 512, 0);    // aspect-preserving
```

Image differences can be computed using `image_difference(a,b)`. This function
performs a simple per-pixel difference and returns the image so that one can
use any metric to determine whether images where different within some threshold.
Optionally, a difference image is returned that highlights the diff.

```cpp
auto a = make_image(...), b = make_image(...); // init images
auto diff = image_difference(a, b);            // image difference
auto display = image_difference(a,nb, true);   // diff display
```

## Procedural images

Yocto/Image defines several procedural images used for both testing and to
quickly create textures for procedural scenes. Testing patterns take as input
the desired image size, the pattern scale, and optionally two colors to
interpolate between. Use `make_grid(...)` for a grid image,
`make_checker(...)` for a checker, `make_bumps(...)` for a bumpy test,
`make_ramp(...)` for a linear ramp, `make_gammaramp(...)` for a ramp with
different gamma settings, `make_uvramp(...)` and `make_uvgrid(...)`
for test images for texture coordinates as either ramps or grids,
`make_colormapramp(...)` for a ramp to test different color maps, and
`make_blackbodyramp(...)` for a ramp with different blackbody temperatures.

Perlin noise patterns can be generated with `make_noisemap(...)`,
`make_fbmmap(...)`, `make_fbmmap(...)` and `make_fbmmap(...)`.
The latter three functions take as input the set of params that control
fractal variations. See [Yocto/Noise](yocto_noise.md) for a description.

```cpp
auto w = 1024, h = 1024;                           // image size
auto scale = float{1};                             // pattern scale
auto c0 = vec4f{0,0,0,1}, c1 = vec4f{0,0,0,1};     // colors
auto i1 = make_grid(w, h, scale, c0, c1);          // grid image
auto i2 = make_checker(w, h, scale, c0, c1);       // checker image
auto i3 = make_bumps(w, h, scale, c0, c1);         // bumps image
auto i4 = make_ramp(w, h, scale, c0, c1);          // ramp image
auto i5 = make_gammaramp(w, h, scale, c0, c1);     // gamma ramp image
auto i6 = make_uvramp(w, h, scale);                // uv ramp image
auto i7 = make_uvgrid(w, h, scale);                // uv grid image
auto i8 = make_colormapramp(w, h, scale);          // color map image

auto t0 = 1000, t1 = 12000;                        // blackbody temperatures
auto b2 = make_blackbodyramp(w, h, scale, t0, t1); // blackbody ramp image

auto noise = vec4f{2, 0.5, 8, 1};                    // noise params
auto n1 = make_noisemap(w, h, scale, noise, c0, c1); // noise image
auto n2 = make_fbmmap(w, h, scale, noise, c0, c1);   // fbm image
auto n3 = make_turbulencemap(w, h, scale, noise, c0, c1);// turbulence image
auto n4 = make_ridgemap(w, h, scale, noise, c0, c1); // ridge image
```

Procedurals skies are generated with
`make_sunsky(img, w, h, elevation, turbidity, sun)`.
The function returns a procedural HDR sky.
The sun position is controlled by its `elevation` that is an angle in `[0,pi/2\]`.
The sky turbidity is controlled by the `turbidity` parameter that is defined in
the range `[1.7,10]`. The `sun` flag determines whether the sun disk is present
in the image. The function support optional parameters to control sun size and
intensity and ground albedo, mostly used for artistic effects.

```cpp
auto sky = make_sunsky(1024, 512, pi/2, 3);        // clear sky
auto sun = make_sunsky(1024, 512, pi/2, 3, true);  // clear sky with sun
auto tur = make_sunsky(1024, 512, pi/2, 10);       // sky with turbidity
```

Use `bump_to_normal(bumps)` to convert a bump map to a normal map, with both
images stored in a linear color spaces.

```cpp
auto bumps = make_bumps(512, 512);    // procedural bump map
auto normal = bump_to_normal(bumps);  // convert bump map to normal map
```

Finally, borders can be added to images using `add_border(img,width,color)`.

```cpp
auto bordered = add_border(img, 1, {0, 0, 0, 1}); // add a thin black border
```

## Low-level operations

Yocto/Image supports versions of the most of the above functions that work
directly on pixel arrays, rather than the image structure. This low-level
interface may be helpful when building applications that mhave their own
image data structure.

In this interface we support arrays of `vec4f` pixels together with arrays of
`vec4b` pixels for 8bit images, the latter with a much smaller number of
functions. In this case, the color space is handled in the application logic.
Conversion between color types is handled with `float_to_byte(pixels)`,
`byte_tyo_float(pixels)`, `rgb_to_srgb(pixels)` and `srgb_to_rgb(pixels)`.

Since most of these functions have signatures similar to the ones above, we
do not document here explicitly.
