# Yocto/ImageIO: Image serialization

Yocto/ImageIO supports loading and saving images from Png, Jpg, Tga, Bmp,
Pnm, Hdr, Exr, Pfm.
Yocto/ImageIO is implemented in `yocto_imageio.h` and `yocto_imageio.cpp`
and depends on `stb_image.h`, `stb_image_write.h`, and `tinyexr.h`.

## Image serialization

Yocto/ImageIO supports reading and writing images. 
Use `load_image(filename, image)` or `image = load_image(filename)` 
to load images, and `save_image(filename, text)` to save it.
Images are stored as `image_data` structs defined in [Yocto/Image](yocto_image.md).
Upon errors, an `io_error` is thrown from all IO functions.
See [Yocto/CommonIO](yocto_commonio.md) for discussion on error handling 
and use without exceptions.

```cpp
auto image = image_data{};
load_image("input_file.png",  image);        // load image
save_image("output_file.png", image);        // save image
auto image1 = load_image("input_file.png");  // alternative load
```
