# Yocto/Image: Image utilities

Yocto/Image is a collection of image utilities useful when writing rendering
algorithms. These include a simple image data structure, color conversion
utilities and tone mapping, loading and saving functionality, and image resizing.
Yocto/Image is implemented in `yocto_image.h` and `yocto_image.cpp`, and
depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h`,
`tinyexr.h` for the image serialization.

## Image representation

Images are stored in a generic container named `image<T>`, where pixels can
be of any type. The convention is Yocto/GL is to use `vec3f` and `vec4f` for
RGB amd RGBA images respectively, or `vec3b` and `vec4b` for 8-bit storage.
The container has a design that is similar in spirit to `std::vector<T>`.
Images are intended to be used as value types. Im memory, image pixels are
stored in 1D array in row-major order.

Images are default-initialized to an empty container, or can be initialized
with a size and either a constant value or a pointer to image data.
Use `img.size()` to get the image width and height as a `vec2i`,
`img.count()` to get the number of pixels, and `img.empty()` to check
whether the image is empty. Images can be resized with `img.resize(size)`,
re-initialized with `img.assign(size, value)` and cleared with `img.clear()`.

```cpp
auto img = image<vec4f>{{512,512}, {1,0,0,1}}; // creates a 512x512 red image
auto size = img.size();             // size == {512,512}
asset(!img.empty());                // check if empty
img.resize({1024,512});             // resize to {1024,512}
img.assign({1024,512}, {0,1,0,1});  // re-initialize as green
```

Image pixels are accessed via the bracket operator using image coordinates
represented as `vec2i`, i.e. `img[{i,j}]`. Pixels can also be accessed with
a sequential index as `img[idx]`. Pixels can iterated with a foreach or
accessed directly by getting a pointer to the underlying data buffer
with `img.data()`.

```cpp
auto img = image<vec4f>{{512,512}, {1,0,0,1}}; // creates a 512x512
img[{10,10}] = {0,1,0,1};     // set pixel at 10,10
for(auto& p : img) print(p);  // iterate over pixels
print(img[1000]);             // sequential access
auto pixels = img.data();     // get pixel buffer
```

## Image sampling

Images are sampled with `eval_image(img,uv)`. Yocto/Image supports sampling
of 3 or 4 channel images with float or byte channels. All results are converted
in float values for increased precision. For byte images, sampling can be
performed in sRGB or linear color space. By default, images are sampled with
bilinear interpolation and tiling, with nearest-neighbor interpolation
and edge-clamp behavior available as overrides.

```cpp
auto hdr = image<vec4f>{...};           // creates a float image
auto ch0 = eval_image(hdr, {0.5,0.5});  // samples the image center
auto ch1 = eval_image(hdr, {0.5,0.5}, false);  // nearest-neighbor interpolation
auto ch2 = eval_image(hdr, {0.5,0.5}, false, false);  // clamp to edge
auto ldr = image<vec4b>{...};           // creates an 8bit image
auto ch0 = eval_image(ldr, {0.5,0.5});  // samples in sRGB, returns linear RGB
auto ch1 = eval_image(ldr, {0.5,0.5}, true);  // treats 8bit values as linear
```

## Color conversions

For convenience, Yocto/Image define functions that convert images between
float and byte channels, and between linear RGB and sRGB representations,
for images that have one, three and four channels. Use `byte_to_float(img)`
and `float_to_byte(img)` for conversion between float and byte images.
Use `rgb_to_srgb(img)` and `srgb_to_rgb(img)` for conversions between
linear and sRGB color spaces.

```cpp
auto img_f = image<vec4f>{...};     // initialize a linear, float image
auto img_b = float_to_byte(img_f);  // convert to 8bit channels
auto img_s = rgb_to_srgb(img);      // convert to sRGB
```

## Tone mapping

HDR images can be tone mapped using `tonemap_image(hdr,e)` that applies
an exposure correction followed by an optional filmic tone curve.
Use `tonemap_image_mt(...)` to quickly tone map large images using
multiple threads.

```cpp
auto hdr = image<vec4f>{...};            // initialize am HDR image
auto ldr = tonemap_imahe(hdr, 0);        // tone mapping
auto flm = tonemap_imahe(hdr, 0, true);  // filmic tone mapping
```

## Color grading

Images can be color graded by applying a set of minimal color grading tools
using `colorgrade_image(img,linear,params)`, in manner similar to
[Hable 2017](http://filmicworlds.com/blog/minimal-color-grading-tools/).
Color grading corrections can be applied on images that are either linear HDR
or non-linear LDR, i.e. sRGB encoded. The results is always an LDR image encoded
in sRGB. Use `colorgrade_image_mt(...)` to quickly tone map large images using
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
auto hdr = image<vec4f>{...};            // initialize am HDR image
auto params = colorgrade_params{};       // default color grading params
params.exposure = 1;                     // set desired params
params.logcontrast = 0.75;               // set desired params
params.tint = compute_white_balance(hdr);// apply white balance
auto ldr = colorgrade_image(hdr, true, params); // color grading
```

## Image resizing

Images are can resized with `resize_image(img,size)`. Just like all other
functions, images are not resized in placed, but a new image is created.
Resizing works for both linear and 8bit images. The size parameter can be
used to perform aspect-preserving resizing by leaving one dimension as zero.

```cpp
auto img = image<vec4f>{...};              // initialize an image
auto res = resize_image(img, {512, 512});  // resizing to fixed size
auto asp = resize_image(img, {512, 0});    // aspect-preserving
```

## Image diffing

Image difference can be computed using `image_difference(a,b)`. This function
performs a simple per-pixel difference and returns the image so that one can
use any metric to determine whether images where different within some threshold.
Optionally, a difference image is returned that highlights the diff.

```cpp
auto a = image<vec4f>{...}, b = image<vec4f>{...}; // init images
auto diff = image_difference(a,b);                 // image difference
auto display = image_difference(a,b, true);        // diff display
```

## Image loading and saving

Images are loaded with `load_image(filename,img,error)` and saved with
`save_image(filename, img, error)`. Both loading and saving take a filename,
an image buffer and return whether or not the image was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
Yocto/Images supports loading and saving to JPG, PNG, TGA, BMP, HDR, EXR.
Yocto/Image supports loading and saving of images that have 1-4 channels and
have either float or byte channels. The number of channel and type is selected
using overloads.

When loading images, the desired number of channels is determined by the
function overload. If the image content does not contain the same number of
channels, these are automatically converted to the desired color format.
Float images are returned encoded in linear RGB color space, while 8bit
images are returned encoded in sRGB. Since color space information is not
present in many images and since it is customary in 3D graphics to misuse
color encoding, e.g. normal maps stored in sRGB, during loading 8bit file
formats are assume to be encoded sRGB, while float formats are assume to
be encoded in linear RGB. When channel type differ between the desired
output and the file one, the appropriate linear-to-sRGB conversion is applied.
Use `is_hdr_filename(filename)` to determine whether a supported format
contains linear float data.

When saving images, float images are expected to be encoded in linear RGB and
are saved as is to file formats capable of storing floats, or converted to
sRGB for file formats that support only that. Byte images are expected
to be encoded is sRGB and saved as is to 8bit file formats or converted to
linear RGB for saving to floating point file formats.

```cpp
auto error = string{};
auto img4f = image<vec4f>{};            // define desired image format
if(!load_image(filename, img4f, error)) // load image as 4-channel linear
  print_error(error);                   // check and print error
auto img4b = image<vec4b>{};            // define desired image format
if(!load_image(filename, img4b, error)) // load image as 4-channel sRGB
  print_error(error);                   // check and print error
auto is_hdr = is_hdr_filename(filename);// check if hdr format
if(is_hdr) ok = load_image(filename, img4f, error); // load as linear
else       ok = load_image(filename, img4b, error); // or sRGB
if(!save_image(filename, img4f, error)) // save image as 4-channel linear
  print_error(error);                   // check and print error
if(!save_image(filename, img4f, error)) // save image as 4-channel linear
  print_error(error);                   // check and print error
```

## Procedural images

Yocto/Image defines several procedural images used for both testing and to
quickly create textures for procedural scenes. Testing patterns take as input
the desired image size, the pattern scale, and optionally two colors to
interpolate between. Use `make_grid(...)` for a grid image,
`make_checker(...)` for a checker, `make_bumps(...)` for a bumpy test,
`make_ramp(...)` for a linear ramp, `make_gammaramp(...)` for a ramp with
different gamma settings, `make_uvramp(...)` and `make_uvgrid(...)`
for test images for texture coordinates as either ramps or grids, and
`make_blackbodyramp(...)` for a ramp with different blackbody temperatures.

Perlin noise patterns can be generated with `make_noisemap(...)`,
`make_fbmmap(...)`, `make_fbmmap(...)` and `make_fbmmap(...)`.
The latter three functions take as input the set of params that control
fractal variations. See [Yocto/Noise](yocto_noise.md) for a description.

```cpp
auto size = vec2i{512, 512};                    // image size
auto scale = float{1};                          // pattern scale
auto c0 = vec4f{0,0,0,1}, c1 = vec4f{0,0,0,1};  // colors
auto img = image<vec4f>{};                      // image buffer
make_grid(img, size, scale, c0, c1);            // grid image
make_checker(img, size, scale, c0, c1);         // checker image
make_bumps(img, size, scale, c0, c1);           // bumps image
make_ramp(img, size, scale, c0, c1);            // ramp image
make_gammaramp(img, size, scale, c0, c1);       // gamma ramp image
make_uvramp(img, size, scale);                  // uv ramp image
make_uvgrid(img, size, scale);                  // uv grid image

auto t0 = 1000, t1 = 12000;                     // blackbody temperatures
make_blackbodyramp(img, size, scale, t0, t1);   // blackbody ramp image

auto noise = vec4f{2, 0.5, 8, 1};               // noise params
make_noisemap(img, size, scale, noise, c0, c1);        // noise image
make_fbmmap(img, size, scale, noise, c0, c1);          // fbm image
make_turbulencemap(img, size, scale, noise, c0, c1);   // turbulence image
make_ridgemap(img, size, scale, noise, c0, c1);        // ridge image
```

Procedurals skies are generated with `make_sunsky(img,size,elevation,turbidity,sun)`.
The function returns an HDR sky generated with the ... algorithm.
The sun position is controlled by its `elevation` that is an angle in `[0,pi/2\]`.
The sky turbidity is controlled by the `turbidity` parameter that is defined in
the range `[1.7,10]`. The `sun` flag determines whether the sun disk is present
in the image. The function support optional parameters to control sun size and
intensity and ground albedo, mostly used for artistic effects.

```cpp
auto img = image<vec4f>{};                         // image buffer
auto sky = make_sunsky(img, {1024,512}, pi/2, 3);  // clear sky
auto sun = make_sunsky(img, {1024,512}, pi/2, 3, true);  // clear sky with sun
auto sky = make_sunsky(img, {1024,512}, pi/2, 10);  // sky with turbidity
```

Use `bump_to_normal(bumps)` to convert a bump map to a normal map, with both
images stored in a linear color spaces.

```cpp
auto bumps = image<vec4f>{};          // image buffer
make_bumps(bumps, {512,512});         // procedural bump map
auto normal = bump_to_normal(bumps);  // convert bump map to normal map
```

Finally, borders can be added to images using `add_border(img,width,color)`.

```cpp
auto img = image<vec4f>{...};                     // image
auto bordered = add_border(img, 1, {0, 0, 0, 1}); // add a thin black border
```
