# Yocto/SceneIO: Scene serialization

Yocto/SceneIO is a collection of functions to load and save images, shapes, and
scene elements, together with path manipulation utilities. The list of 
supported file formats in included in the sections below.
Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`,
and depends on `json.hpp` for Json serialization, and `stb_image.h`,
`stb_image_write.h` and `tinyexr.h` for the image serialization, and
`cgltf.h` and `cgltf_write.h` for glTF support.

## Image serialization

Use `ok = load_image(filename, image, error)` to load images 
and `ok = save_pbrt(filename, image, error)` to save them.  
Both loading and saving take a filename, a reference to an image, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
The library supports loading and saving to JPG, PNG, TGA, BMP, HDR, EXR.

```cpp
auto error = string{};
auto image = image_data{};              // define an image
if(!load_image(filename, image, error)) // load image
  handle_error(error);                   // check and print error
if(!save_image(filename, image, error)) // save image
  handle_error(error);                   // check and print error
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

Use `ok = load_shape(filename, shape, error)` to load shapes 
and `ok = save_shape(filename, shape, error)` to save them.  
Both loading and saving take a filename, a reference to a model, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
The library supports loading and saving to Ply, Obj, Stl.

```cpp
auto shape = shape_data{};                   // shape
auto error = string{};                       // error buffer
if(!load_shape(filename, shape, error))      // load shape
  handle_error(error);
if(!save_shape(filename, shape, error))      // save shape
  handle_error(error);
```

During loading, the library chooses whether to store triangles or quads.
By default, texture coordinates are flipped vertically to match the convention 
of OpenGL texturing; this can be disabled by setting the `flipv` flag.

## Face-Varying shape serialization

Use `ok = load_fvshape(filename, fvshape, error)` to load face-varying shapes 
and `ok = save_fvshape(filename, fvshape, error)` to save them.  
Both loading and saving take a filename, a reference to a shape, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
The library supports loading and saving to Ply, Obj, Stl.

```cpp
auto fvshape = fvshape_data{};               // face-varying shape
auto error = string{};                       // error buffer
if(!load_fvshape(filename, fvshape, error))  // load shape
  handle_error(error);
if(!save_fvshape(filename, fvshape, error))  // save shape
  handle_error(error);
```

Face-varying shapes are coerced to quads, regardless of the polygon type stored 
in the file. By default, texture coordinates are flipped vertically to match the 
convention of OpenGL texturing; this can be disabled by setting the `flipv` flag.

## Scene serialization

Use `ok = load_scene(filename, scene, error)` to load scenes 
and `ok = save_scene(filename, scene, error)` to save them.  
Both loading and saving take a filename, a reference to a scene, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
The library supports loading and saving to a custom Json format,
Obj, glTF, Pbrt, and all shape file formats.

```cpp
auto scene = scene_data{};               // scene
auto error = string{};                   // error buffer
if(!load_scene(filename, scene, error))  // load scene
  handle_error(error);
if(!save_scene(filename, scene, error))  // save scene
  handle_error(error);
```

## Texture serialization

Use `ok = load_texture(filename, texture, error)` to load textures 
and `ok = save_texture(filename, texture, error)` to save them.  
Both loading and saving take a filename, a reference to a texture, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
Yocto/SceneIO supports loading and saving to JPG, PNG, TGA, BMP, HDR, EXR.

```cpp
auto error = string{};
auto texture = texture_data{};              // define a texture
if(!load_texture(filename, texture, error)) // load texture
  handle_error(error);                       // check and print error
if(!save_texture(filename, texture, error)) // save texture
  handle_error(error);                       // check and print error
```

When loading textures, the color space and bit depth is detected by the filename.
Texture are loaded as float for float file formats, and as bytes for 8bit file
formats. Since color space information is not present in many images,
we assume that byte textures are encoded in sRGB, and float texture in linear RGB.
Use `is_hdr_filename(filename)` to determine whether a supported format
contains linear float data.
When saving images, pixel values are converted to the color space supported
by the chosen file format.

## Text and binary serialization

Use `ok = load_text(filename, text, error)` to load text files 
and `ok = save_text(filename, text, error)` to save them.  
Use `ok = load_binary(filename, data, error)` to load binary files 
and `ok = save_binary(filename, data, error)` to save them. 
Both loading and saving take a filename, a reference to a string or data buffer, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.
Text is stored as a string and binary data is stored as an array of bytes.
Use without exception is supported as described above.

```cpp
auto text = string{};
auto error = string{};
if(!load_text("input_file.txt",  text, error))          // load text
  handle_error(error);
if(!save_text("output_file.txt", text, error))          // save text
  handle_error(error);
auto data = vector<byte>{};
load_binary("input_file.bin",  data);        // load data
  handle_error(error);
save_binary("output_file.bin", data);        // save data
  handle_error(error);
```

## Path manipulation utilities

Yocto/CommonIO contains several helper function to manipulate paths. Paths
are encoded in UTF8 across the library and these functions make it easier to
handle UTF8-encoded paths across operating systems, by wrapping 
`std::filesystem` with a string interface.

Use `path_dirname(filename)`, `path_extension(filename)`,  
`path_filename(filename)`, `path_basename(fillename)`
to extract the directory, extension, filename and basename from a path.
Use `path_join(patha,pathb)` to joins paths and
`replace_extension(filename,ext)` to replace a path extension.
Use `path_exists(filename)` to check if a path exists and
`path_isdir(filename)` and `path_isfile(filename)` to check whether
it is a directory ot file respectively.
Use `ok = list_directory(dirname, error)` to list directory contents, 
`ok = make_directory(dirname, error)` to create a directory, and
`path_current()` to get the current directory.
