# Yocto/SceneIO: Scene serialization

Yocto/SceneIO supports loading and saving scenes Yocto/Scene models from
Ply, Obj, Pbrt, glTF and a custom Json format.
Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`, and depends on `cgltf.h`.

## Serialization formats

Yocto/SceneIO supports loading and saving to Ply, Obj, Pbrt, glTF,
and a custom Json format. For the standard formats, loading is best effort,
since scene data is transformed from the formats' scene models to the
Yocto/Scene model.

The custom Json format is a serialization of the internal properties for
most scene objects, with a few conventions taken for extensibility.
Scene's object arrays are represented as dictionaries in Json with the
objects' names used as keys. This ensure proper reference semantic and
allows for more extensibility in the future, but it also means that
object order is not preserved during serialization.

Cameras, materials, instances and environments are represented directly
in the Json scene, while textures and shapes are serialized using
standard image and geometry formats. By convention, scenes are
stored as a single Json format for the scene structure. Textures are
stored in the `textures` directory with the name of the texture as filename,
while the extension is determined by checking th available files.
Shapes are stored in the `shapes` directory with name of the shape as filename,
while the extension is determined by checking th available files.

## Loading and saving scenes

Scenes are loaded with `load_scene(filename, scene, error, progress)` and
saved with `save_scene(filename, scene, error, progress)`.
Both loading and saving take a filename, a scene pointer and return
whether or not the scene was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
The functions take a progress callback as an optional parameter,
that is called as scene loading progresses.

```cpp
auto scene = new scene_model{};                    // scene
auto progress = [](const string& message,          // progress callback
                   int current, int total) {
  print_info(message, current, total);
};
auto error = string{};                             // error buffer
if(!load_scene(filename, scene, error, progress))  // load scene
  print_error(error);
if(!save_scene(filename, scene, error, progress))  // save scene
  print_error(error);
```
