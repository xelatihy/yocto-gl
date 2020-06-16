# Yocto/Image: Image utilities

Yocto/Image is a collection of image utilities useful when writing rendering
algorithms. These include a simple image data structure, color conversion
utilities and tone mapping. We provide loading and saving functionality for
images and support PNG, JPG, TGA, BMP, HDR, EXR formats.
Yocto/Image is implemented in `yocto_image.h` and `yocto_image.cpp`, and
depends on stb_image.h, stb_image_write.h, stb_image_resize.h,
tinyexr.h for the image serialization

## Images

Yocto/Math contains a simple image container that can be used to store
generic images. The container is similar in spirit to `std::vector`.
We provide only minimal image functions including lookup and sampling.

## Image Utilities

Yocto/Image supports a very small set of color and image utilities including
color utilities, example image creation, tone mapping, image resizing, and
sunsky procedural images. Yocto/Image is written to support the need of a
global illumination renderer, rather than the need of generic image editing.
We support 4-channels float images (assumed to be in linear color) and
4-channels byte images (assumed to be in sRGB).

1. store images using the image<T> structure
2. load and save images with `load_image()` and `save_image()`
3. resize images with `resize()`
4. tonemap images with `tonemap()` that convert from linear HDR to
   sRGB LDR with exposure and an optional filmic curve
5. make various image examples with the `make_proc_image()` functions
6. create procedural sun-sky images with `make_sunsky()`
7. many color conversion functions are available in the code below
