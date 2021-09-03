# Yocto/SceneIO: Scene serialization

Yocto/SceneIO is a collection of functions to load and save images, shapes, and
scene elements, together with path manipulation utilities.
Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`,
and depends on `json.hpp` for Json serialization, and `stb_image.h`,
`stb_image_write.h` and `tinyexr.h` for the image serialization.

## Errors handling in IO functions

IO functions in Yocto/GL have a dual interface to support their use with and 
without exceptions. IO functions with exception are written as 
`load_<type>(filename, data, <options>)` and 
`save_<type>(filename, data, <options>)` where `<type>` is the data type 
read or written, and `<options>` is an optional list of IO options. 
A convenience shortcut is provided for all loading functions that returns the 
data directly, written as `data = load_<type>(filename, <options>)`.
Upon errors, an `io_error` is thrown from all IO functions.
This makes the error-handling code more uniform across the library. 
`io_error` has fields for `filename` and `message` that can retrieved directly, 
and a message for users that can be retrieved with the `.what()` method.

```cpp
auto text = string{};
try {
  load_text("input_file.txt",  text);   // load text
  text = load_text("input_file.txt");   // alternative load
  save_text("output_file.txt", text);   // save text
} catch(const io_error& error) {
  print_info(error.what()); exit(1);    // handle error
}
```

IO functions without exceptions are written as
`load_<type>(filename, data, error, <options>)` and 
`save_<type>(filename, data, error, <options>)` where `error` is a string 
that contains a message if an error occurred. These functions return a boolean 
flag that indicates whether the operation succeeded.  

```cpp
auto text = string{};
auto error = string{};
if(!load_text("input_file.txt",  text))      // load text
  print_info(error); exit(1);                // handle error
if(!save_text("output_file.txt", text))      // save text
  print_info(error); exit(1);                // handle error
```

## Image serialization

Images are loaded with `load_image(filename, image)` or 
`image = load_image(filename)` and saved with `save_image(filename, text)`.
Images are stored as `image_data` structs defined in [Yocto/Image](yocto_image.md).
Upon errors, an `io_error` is thrown from all IO functions.
Yocto/SceneIO supports loading and saving to JPG, PNG, TGA, BMP, HDR, EXR.

```cpp
auto image = image_data{};
load_image("input_file.png",  image);        // load image
save_image("output_file.png", image);        // save image
auto image1 = load_image("input_file.png");  // alternative load
```

When using without exceptions, images are loaded with 
`load_image(filename, img, error)` and saved with
`save_image(filename, img, error)`. Both loading and saving take a filename,
an image buffer and return whether or not the image was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.

```cpp
auto error = string{};
auto image = image_data{};              // define an image
if(!load_image(filename, image, error)) // load image
  print_error(error);                   // check and print error
if(!save_image(filename, image, error)) // save image
  print_error(error);                   // check and print error
```

When loading images, the color space is detected by the filename and stored
in the returned image. Since color space information is not
present in many images and since it is customary in 3D graphics to misuse
color encoding, e.g. normal maps stored in sRGB, during loading 8bit file
formats are assume to be encoded sRGB, while float formats are assumed to
be encoded in linear RGB.
Use `is_hdr_filename(filename)` to determine whether a supported format
contains linear float data.
When saving images, pixel values are converted to the color space supported
by the chosen file format.

## Shape serialization

Shapes are loaded with `load_shape(filename, shape, error)` and
saved with `save_shape(filename, shape, error)`.
Face-varying shapes with `load_fvshape(filename, shape, error)`
and saved with `save_fvshape(filename, shape, error)`.
Both loading and saving take a filename, a scene reference and return
whether or not the scene was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
Yocto/SceneIO supports loading and saving to Ply, Obj, Stl.

For indexed meshes, the type of elements is determined during loading,
while face-varying coerce all face to quads. By default, texture coordinates
are flipped vertically to match the convention of OpenGL texturing;
this can be disabled by setting the `flipv` flag.

```cpp
auto shape = shape_data{};                   // shape
auto error = string{};                       // error buffer
if(!load_shape(filename, shape, error))      // load shape
  print_error(error);
if(!save_shape(filename, shape, error))      // save shape
  print_error(error);

auto fvshape = fvshape_data{};               // face-varying shape
auto error = string{};                       // error buffer
if(!load_fvshape(filename, fvshape, error))  // load shape
  print_error(error);
if(!save_fvshape(filename, fvshape, error))  // save shape
  print_error(error);
```

## Scene serialization

Scenes are loaded with `load_scene(filename, scene, error, progress)` and
saved with `save_scene(filename, scene, error, progress)`.
Both loading and saving take a filename, a scene reference and return
whether or not the scene was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
Yocto/SceneIO supports loading and saving to a custom Json format,
Obj, glTF, Pbrt, and all shape file formats.

```cpp
auto scene = scene_data{};               // scene
auto error = string{};                   // error buffer
if(!load_scene(filename, scene, error))  // load scene
  print_error(error);
if(!save_scene(filename, scene, error))  // save scene
  print_error(error);
```

## Texture serialization

Textures are loaded with `load_texture(filename, texture, error)` and saved with
`save_texture(filename, texture, error)`. Both loading and saving take a filename,
a texture and return whether or not the image was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
Yocto/SceneIO supports loading and saving to JPG, PNG, TGA, BMP, HDR, EXR.

When loading textures, the color space and bit depth is detected by the filename.
Texture are loaded as float for float file formats, and as bytes for 8bit file
formats. Since color space information is not present in many images,
we assume that byte textures are encoded in sRGB, and float texture in linear RGB.
Use `is_hdr_filename(filename)` to determine whether a supported format
contains linear float data.
When saving images, pixel values are converted to the color space supported
by the chosen file format.

```cpp
auto error = string{};
auto texture = texture_data{};              // define a texture
if(!load_texture(filename, texture, error)) // load texture
  print_error(error);                       // check and print error
if(!save_texture(filename, texture, error)) // save texture
  print_error(error);                       // check and print error
```
