# Yocto/Img

Utilities for loading/saving images, mostly to avoid linking
problems and code duplication across other yocto libraries. Functions as
a wrapper to other libraries.

This library depends in yocto_math.h, stb_image.h, stb_image_write.h
stb_image_resize.h.


## Usage

1. load images with `load_image()`
2. save images with `save_image()`
2. resize images with `resize_image()`


## History

- v 0.3: use reference interface for images
- v 0.2: added whole image functions
- v 0.1: initial release

## Namespace yimg

Utilitied for reading and writing images.

### Function load_image4b()

~~~ .cpp
ym::image4b load_image4b(const std::string& filename);
~~~

Loads an ldr image.

### Function load_image4f()

~~~ .cpp
ym::image4f load_image4f(const std::string& filename);
~~~

Loads an hdr image.

### Function save_image4b()

~~~ .cpp
bool save_image4b(const std::string& filename, const ym::image4b& img);
~~~

Saves an ldr image.

### Function save_image4f()

~~~ .cpp
bool save_image4f(const std::string& filename, const ym::image4f& img);
~~~

Saves an hdr image.

### Function load_image()

~~~ .cpp
byte* load_image(const std::string& filename, int& w, int& h, int& ncomp);
~~~

Loads an ldr image.

### Function load_imagef()

~~~ .cpp
float* load_imagef(const std::string& filename, int& w, int& h, int& ncomp);
~~~

Loads an hdr image.

### Function save_image()

~~~ .cpp
bool save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr);
~~~

Saves an ldr image. Uses extension to determine which format to load.
Suppoted formats are PNG.

### Function save_imagef()

~~~ .cpp
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr);
~~~

Saves an hdr image. Uses extension to determine which format to load.
Suppoted formats are HDR.

### Function load_image_from_memory()

~~~ .cpp
byte* load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp);
~~~

Loads an image from memory.

### Function load_imagef_from_memory()

~~~ .cpp
float* load_imagef_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp);
~~~

Loads an image from memory.

### Function load_image4b_from_memory()

~~~ .cpp
ym::image4b load_image4b_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp);
~~~

Loads an image from memory.

### Function load_image4f_from_memory()

~~~ .cpp
ym::image4f load_image4f_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp);
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

