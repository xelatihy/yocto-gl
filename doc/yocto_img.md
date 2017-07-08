# Yocto/Img

Utilities for loading/saving images, mostly to avoid linking
problems and code duplication across other yocto libraries. Functions as
a wrapper onto other libraries.

This library depends in yocto_math.h, stb_image.h, stb_image_write.h
stb_image_resize.h.


## Usage

1. load images with `load_image()`
2. save images with `save_image()`
2. resize images with `resize_image()`


## History

- v 0.1: initial release

## Namespace yimg

Utilitied for reading and writing images.

### Function load_image()

~~~ .cpp
void load_image(
    const std::string& filename, int& w, int& h, int& ncomp, byte*& ldr);
~~~

Loads an ldr image.

### Function load_image()

~~~ .cpp
void load_image(
    const std::string& filename, int& w, int& h, int& ncomp, float*& hdr);
~~~

Loads an hdr image.

### Function load_image()

~~~ .cpp
void load_image(const std::string& filename, int& w, int& h, int& ncomp,
    float*& hdr, byte*& ldr);
~~~

Loads an image. Uses extension to determine which format to load.
Suppoted formats are PNG, JPG, HDR.

### Function load_image()

~~~ .cpp
void load_image(const std::string& filename, int& w, int& h, int& ncomp,
    std::vector<float>& hdr, std::vector<byte>& ldr);
~~~

Loads an image.

### Function save_image()

~~~ .cpp
void save_image(const std::string& filename, int width, int height, int ncomp,
    const float* hdr);
~~~

Saves an hdr image. Uses extension to determine which format to load.
Suppoted formats are HDR.

### Function save_image()

~~~ .cpp
void save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr);
~~~

Saves an ldr image. Uses extension to determine which format to load.
Suppoted formats are PNG.

### Function load_image_from_memory()

~~~ .cpp
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp, float*& hdr, byte*& ldr);
~~~

Loads an image from memory.

### Function load_image_from_memory()

~~~ .cpp
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp, std::vector<float>& hdr,
    std::vector<byte>& ldr);
~~~

Loads an image from memory.

### Function resize_image()

~~~ .cpp
void resize_image(int width, int height, int ncomp, const float* img,
    int res_width, int res_height, float* res_img);
~~~

Resize image.

### Function resize_image()

~~~ .cpp
void resize_image(int width, int height, int ncomp, const byte* img,
    int res_width, int res_height, byte* res_img);
~~~

Resize image.

### Function resize_image()

~~~ .cpp
void resize_image(int width, int height, int ncomp,
    const std::vector<float>& img, int res_width, int res_height,
    std::vector<float>& res_img);
~~~

Resize image.

### Function resize_image()

~~~ .cpp
void resize_image(int width, int height, int ncomp,
    const std::vector<byte>& img, int res_width, int res_height,
    std::vector<byte>& res_img);
~~~

Resize image.

