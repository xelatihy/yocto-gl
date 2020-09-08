# Yocto/ModelIO: Serialization for Obj, Ply and Pbrt models

Yocto/ModelIO is a collection of utilities for loading and saving scenes
and meshes in Ply, Obj and Pbrt formats.
Yocto/ModelIO is implemented in `yocto_modelio.h` and `yocto_modelio.cpp`.

## Ply models

The Ply file format is a generic file format used to serialize meshes and point
clouds. To use this library is helpful to understand the basic of the Ply
file format for example from the
[Ply Wikipedia page](<https://en.wikipedia.org/wiki/PLY_(file_format)>).

Yocto/ModelIO represents Ply data with the `ply_model` struct.
The internal representation matches the structure of a Ply file and
can be accessed directly if desired. The Ply model is defined as an array
of Ply elements, which in turn are defined as arrays of Ply properties.
Ply properties can contain most C data types.
All elements and properties are owned by the main `ply_model`.
Yocto/ModelIO provides several functions to read and write Ply data
whose use is preferred over direct data access.

Use `load_ply(filename, ply, error)` to load Ply files and
`save_ply(filename, ply, error)` to save them.  
Both loading and saving take a filename, a pointer to a Ply model,
and returns whether or not the file was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.

```cpp
auto ply = new ply_model{};             // ply model buffer
auto error = string{};                  // error buffer
if(!load_ply(filename, ply, error))     // load ply
  print_error(error);                   // check and print error
if(!save_ply(filename, ply, error))     // save ply
  print_error(error);                   // check and print error
```

## Ply reading

Yocto/ModelIO defines several functions to make it easy to extract data
from Ply files. Since Ply properties cane be stored with many different
C types, the convenience functions convert the various underlying
representations to the requested one.

Use `has_property(ply,element,property)` to check whether the model has
a property named `property` in an element named `element`. Use
`get_property(ply,element,property)` to get that property.

Use `get_value(ply, element, property, values)` to get the property values
and stored them in the array `values`. The function returns whether or not
the property was present. Since Ply properties are often stored together,
e.g. xyz coordinates, Yocto/ModelIO supports queries that take arrays of
properties and returns values packed in `vecXf`.
Use `get_values(ply, element, properties, values)` to read arrays of
properties.

For list properties, Yocto/ModelIO supports reading properties as
arrays of arrays of dynamic size with `get_lists(ply,element,property,lists)`.
A faster, but harder to use, method is to get lists sizes and values as
separate arrays, where list values are packed together to avoid small memory
allocations. Use `get_list_sizes(ply, element, property, sizes)` for sizes
and `get_list_values(ply, element, property, values)` for values.

```cpp
auto ply = new ply_model{};             // ply model buffer
auto error = string{};                  // error buffer
load_ply(filename, ply, error);         // load ply

auto radius = vector<float>{};                   // property buffer
if(!get_value(ply, "vertex", "radius", radius))  // read property
  print_error("missing radius");                 // error if missing
auto positions = vector<vec3f>{};                // properties buffer
if(!get_values(ply, "vertex", {"x","y","z"}, positions)) // read properties
  print_error("missing positions");              // error if missing

auto faces = vector<vector<int>>{};              // list property buffer
if(!get_lists(ply, "face", "indices", faces))    // read lists
  print_error("missing faces");                  // error if missing

auto faces_sizes = vector<int>{};                // list property sizes
auto faces_values = vector<int>{};               // list property values
if(!get_list_sizes(ply, "face", "indices", faces_sizes)) // read lists sizes
  print_error("missing faces");                  // error if missing
if(!get_list_values(ply, "face", "indices", faces_values)) // read lists values
  print_error("missing faces");                  // error if missing
```

Yocto/Shape defines convenience functions to read the most used properties
of meshes, using standard element and property names.
For vertex properties, use `get_positions(ply, positions)`,
`get_normals(ply, normals)`, `get_texcoords(ply, texcoords, flipv)`,
`get_colors(ply, colors)`, and `get_radius(ply, radius)` to read positions,
normals, texcoords, colors and radius if present. Texture coordinates can be
optionally flipped vertically. For shape elements, use
`get_points(ply, points)`, `get_lines(ply, lines)`,
`get_triangles(ply, triangles)`, and `get_quads(ply, quads)` , to read
points, lines, triangles and quads. Note that since Ply support arbitrary
polygons and polylines, these functions tesselate the Ply polygons into
the desired element type, for now using a simple fan-like algorithm that
works only for convex elements. Use `has_quads(ply)` to check whether
the Ply data has quads.

```cpp
auto ply = new ply_model{};             // ply model buffer
auto error = string{};                  // error buffer
load_ply(filename, ply, error);         // load ply

auto positions = vector<vec3f>{};       // vertex properties buffers
auto normals   = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto colors    = vector<vec3f>{};
auto radius    = vector<float>{};
get_positions(ply, positions);          // read vertex props.
get_normals(ply, normals);
get_texcoords(ply, texcoords, false);   // last params, flips y is desired
get_colors(ply, colors);
get_radius(ply, radius);

auto points    = vector<vec4i>{};       // shape elements buffers
auto lines     = vector<vec4i>{};
auto triangles = vector<vec4i>{};
auto quads     = vector<vec4i>{};
get_points(ply, points);
get_lines(ply, lines);
if(has_quads(ply)) get_quads(ply, quads);
else get_triangles(ply, triangles);
```

## Ply writing

Yocto/ModelIO defines several functions to make it easy to fill Ply data
to save to file. Since Ply properties cane be stored with many different
C types, the convenience functions maintain the same data type of the
data passed in without performing any conversion.

Use `add_value(ply, element, property, values)` to add the property values
stored them in the array `values`. The function returns whether or not
the property was successfully added. Since Ply properties are often stored together,
e.g. xyz coordinates, Yocto/ModelIO supports adding arrays of properties whose
values are packed in `vecXf`.
Use `add_values(ply, element, properties, values)` to add arrays of
properties.

For list properties, Yocto/ModelIO supports adding list properties as
arrays of arrays of dynamic size with `add_lists(ply,element,property,lists)`.
A faster, but harder to use, method is to add lists sizes and values as
separate arrays, where list values are packed together to avoid small memory
allocations. Use `add_lists(ply, element, property, sizes, values)`.
Finally, Yocto/ModelIO supports adding lists of fixed lengths, where the
parameters are packed into `vecXi`,
with `add_lists(ply,element,property,values)`.

```cpp
auto ply = new ply_model{};             // ply model buffer

auto radius = vector<float>{...};                // property buffer
if(!add_value(ply, "vertex", "radius", radius))  // add property
  print_error("error in radius");                // error if missing
auto positions = vector<vec3f>{...};             // properties buffer
if(!add_values(ply, "vertex", {"x","y","z"}, positions)) // add properties
  print_error("missing positions");              // error if missing

auto faces = vector<vector<int>>{...};           // list property buffer
if(!add_lists(ply, "face", "indices", faces))    // add lists
  print_error("missing faces");                  // error if missing

auto faces_sizes = vector<int>{...};             // list property sizes
auto faces_values = vector<int>{...};            // list property values
if(!add_lists(ply, "face", "indices", faces_sizes, face_values)) // add lists
  print_error("missing faces");                  // error if missing

auto triangles = vector<vec3i>{...};             // fixed length list property
if(!add_lists(ply, "face", "indices", triangles))// add lists
  print_error("missing faces");                  // error if missing

auto error = string{};                  // error buffer
save_ply(filename, ply, error);         // save ply
```

Yocto/ModelIO defines convenience functions to add the most used properties
of meshes, using standard element and property names.
For vertex properties, use `add_positions(ply, positions)`,
`add_normals(ply, normals)`, `add_texcoords(ply, texcoords, flipv)`,
`add_colors(ply, colors)`, and `add_radius(ply, radius)` to read positions,
normals, texcoords, colors and radius if present. Texture coordinates can be
optionally flipped vertically. For shape elements, use
`add_points(ply, points)`, `add_lines(ply, lines)`,
`add_triangles(ply, triangles)`, and `add_quads(ply, quads)`, to add
points, lines, triangles and quads. Use `add_faces(ply, faces)` to add
arbitrary polygonal faces.

```cpp
auto ply = new ply_model{};             // ply model buffer

auto positions = vector<vec3f>{...};    // vertex properties buffers
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>{...};
auto colors    = vector<vec3f>{...};
auto radius    = vector<float>{...};
add_positions(ply, positions);          // add vertex props.
add_normals(ply, normals);
add_texcoords(ply, texcoords, false);   // last params, flips y is desired
add_colors(ply, colors);
add_radius(ply, radius);

auto points    = vector<vec4i>{...};    // shape elements buffers
auto lines     = vector<vec4i>{...};
auto triangles = vector<vec4i>{...};
auto quads     = vector<vec4i>{...};
add_points(ply, points);                // add shape elements
add_lines(ply, lines);
add_triangles(ply, triangles);
add_quads(ply, quads);

auto error = string{};                  // error buffer
save_ply(filename, ply, error);         // save ply
```

## Obj models

The Obj file format is a file format used to serialize meshes and materials.
To use this library is helpful to understand the basic of the Obj
file format for example from the
[Obj Wikipedia page](<https://en.wikipedia.org/wiki/OBJ_(file_format)>).
Obj files come in pairs, `.obj` for shapes and `.mtl` for materials.

Yocto/ModelIO represents Obj data with the `obj_scene` struct.
Obj models are defined as collections of shapes and materials.
Obj shapes use a face-varying representation that has vertex positions,
normals and texture coordinates, with their their own topology.
In a shape, each face is tagged with a material used for that face.
Yocto/ModelIO provides direct access to these tagged shapes by inspecting the
`obj_shape` properties.

In Yocto/Obj, materials are represented by the `obj_material` struct.
Each material is a collection of values and textures that specify the material
lobes, like emission, diffuse, specular, reflection, etc. Each value
has a corresponding texture stored by saving its filename and lookup properties
in `obj_texture_info`. The meaning of these parameters often depend on the
application since Obj was used in many different contexts by reinterpreting the
material values.

Yocto/ModelIO defines three extensions to the Obj file format. Materials have
additional parameters that specify a PBR parametrization similar to glTF, which
is common used in most engines today. Yocto/ModelIO adds another file,
namely `.objx`, that stores cameras and environment maps, respectively as
`obj_camera` and `obj_environment`. This addition makes the extended Obj format
capable of storing full scenes.

The Obj model is defined as an array of objects of the types defined above.
Obj objects are pointers owned by the main `obj_scene`.
Objects properties can be read and written directly from the model data,
and are documented in the header file for now.
For shapes, Yocto/ModelIO provides several functions to read and write Obj
shapes, with a simpler interface than accessing data directly.

```cpp
auto obj = new obj_scene{...};             // obj model buffer
for(auto shape : obj->shapes)              // access shapes
  print_info(shape->name);                 // access shape properties
for(auto material : obj->material)         // access materials
  print_info(material->diffuse);           // access material properties
for(auto material : obj->material)         // access materials
  print_info(material->diffuse_tex);       // access material textures
for(auto camera : obj->cameras)            // access cameras [extension]
  print_info(camera->frame);               // access camera properties
for(auto environment : obj->environments)  // access environments [extension]
  print_info(environment->emission);       // access environment properties
```

Use `load_obj(filename, obj, error)` to load Obj files and
`save_obj(filename, obj, error)` to save them.  
Both loading and saving take a filename, a pointer to a Obj model,
and returns whether or not the file was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.

```cpp
auto obj = new obj_scene{};             // obj model buffer
auto error = string{};                  // error buffer
if(!load_obj(filename, obj, error))     // load obj
  print_error(error);                   // check and print error
if(!save_obj(filename, obj, error))     // save obj
  print_error(error);                   // check and print error
```

## Obj reading

Obj is a face-varying format and that geometry representation is maintained
in `obj_shape`. Yocto/ModelIO provides easier accessed to Obj shape data,
both as indexed meshes and as face-varying meshes.

To get data as a standard indexed meshes, use
`get_triangles(obj, triangles, <vertex>, <materials>)`,
`get_quads(obj, quads, <vertex>, <materials>)`,
`get_lines(obj, lines, <vertex>, <materials>)`, and
`get_points(obj, points, <vertex>, <materials>)`,
to read triangles, quads, lines and points respectively.
In these functions, vertex data is comprised of positions, normals and texture
coordinated stored as separate arrays.
Note that in these functions, vertices may end up being duplicated when
going from the face-varying representation to an indexed mesh.
Material data is comprised of a list of material names used in the shape and
the per-element indices to the material arrays.
Since Obj stored faces as polygons, these functions are performing a tesselation
when necessary that for now work only for convex shapes. You can check whether
a shape contains quads with `has_quads(shape)`.

In some cases, it may be desireable to extract the shape elements corresponding
to a single material, for example for use in renderers that support a single
shader per shape. Yocto/ModelIO defines convenience functions for this case,
that are overrides of the previous `get_<element>(...)` functions,
but differ in that they take a material is as input, instead of returning
materials tags.

```cpp
auto obj = new obj_scene{};             // obj model buffer
auto error = string{};                  // error buffer
load_obj(filename, obj, error);         // load obj
auto shape = obj->shapes.front();       // get shape

auto triangles = vector<vec3i>{};       // element data
auto quads     = vector<vec4i>{};
auto positions = vector<vec3f>{};       // vertex properties
auto normals   = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto materials = vector<string>{};      // materials
auto ematerials = vector<int>{};        // per-face material ids
if(has_quads(shape)) {
  get_triangles(shape, triangles,       // read as triangles
    positions, normals, texcoords,
    materials, ematerials, false);      // flip texcoords if desired
} else {
  get_quads(shape, quads,               // read as quads
    positions, normals, texcoords,
    materials, ematerials, false);      // flip texcoords if desired
}

auto material_id = 0;                     // material kd to extract to
if(has_quads(shape)) {
  get_triangles(shape, material_id,       // read as triangles for material id
    triangles, positions,
    normals, texcoords, false);           // flip texcoords if desired
} else {
  get_quads(shape, material_id,           // read as quads for material 0
    quads, positions,
    normals, texcoords, false);           // flip texcoords if desired
}
```

Yocto/ModelIO supports also reading Obj shapes as face-varying quads
with `get_fvquads(...)`.

```cpp
auto obj = new obj_scene{};             // obj model buffer
auto error = string{};                  // error buffer
load_obj(filename, obj, error);         // load obj
auto shape = obj->shapes.front();       // get shape

auto quadspos  = vector<vec4i>{};       // face-varying element data
auto quadsnorm = vector<vec4i>{};
auto quadsuv   = vector<vec4i>{};
auto positions = vector<vec3f>{};       // vertex properties
auto normals   = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto materials = vector<string>{};      // materials
auto ematerials = vector<int>{};        // per-face material ids
get_fvquads(shape,                      // read as face-varying quads
  quadspos, quadsnorm, quadsuv,
  positions, normals, texcoords,
  materials, ematerials, false);        // flip texcoords if desired
```

## Obj writing

To save an Obj, create a scene and add objects to it using
`add_camera(obj)`, `add_material(obj)`, `add_environment(obj)`
and `add_shape(obj)` for camera, materials, environments and shapes
respectively. For all objects, set the objects' properties directly.

For shapes, Yocto/ModelIO defines convenience functions that take either
indexed mesh or face-varying meshes as input and create the appropriate
Obj shape elements.
Use `set_triangles(shape, triangles, <vertex>, <materials>)`,
`set_quads(shape, quads, <vertex>, <materials>)`,
`set_lines(shape, lines, <vertex>, <materials>)`,
`set_points(shape, points, <vertex>, <materials>)` to set shapes as an indexed
mesh of triangles, quads, lines or points respectively.
In these functions, vertex data is comprised of positions, normals and texture
coordinated stored as separate arrays.
Material data is only represented as tags and can be left
empty if only one material is used. To set material names
use `set_materials(shape, materials)`.

```cpp
auto obj = new obj_scene{};             // obj model buffer

auto camera = add_camera(obj);          // add camera
camera->name = "camera";                // set camera name
camera->frame = identity3x4f;           // set camera properties
auto environment = add_environment(obj);// add environment
environment->name = "environment";      // set environment name
environment->emission = {1,1,1};        // set environment properties
auto material = add_material(obj);      // add material
material->name = "camera";              // set material name
material->diffuse = {1,0,0};            // set material properties

auto triangles = vector<vec3i>{...};    // element data
auto positions = vector<vec3f>{...};    // vertex properties
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>{...};

auto shape = add_shape(obj);            // add shape
shape->name = "shape";                  // set shape name
set_triangles(shape, triangles,         // set shape geometry
  positions, normals, texcoords);
set_materials(shape, {material});       // set shape material;

auto error = string{};                  // error buffer
save_obj(filename, obj, error);         // save obj
```

Yocto/ModelIO supports also writing of face-varying shapes with
`set_fvquads(...)` with an API similar to above.

```cpp
auto obj = new obj_scene{};             // obj model buffer

auto quadspos  = vector<vec4i>{...};    // face-varying element data
auto quadsnorm = vector<vec4i>{...};
auto quadsuv   = vector<vec4i>{...};
auto positions = vector<vec3f>{...};    // vertex properties
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>....{};

auto shape = add_shape(obj);            // add shape
shape->name = "shape";                  // set shape name
set_fvquads(shape, quadspos,            // set shape geometry
  quadsnorm, quadstexcoord,
  positions, normals, texcoords);
set_materials(shape, {material});       // set shape material;

auto error = string{};                  // error buffer
save_obj(filename, obj, error);         // save obj
```

## Pbrt models

The Pbrt file format is a scene representation suitable for realistic rendering
and implemented by the Pbrt renderer.
To use this library is helpful to understand the basic of the Pbrt
file format for example from the
[Pbrt format documentatioon](https://www.pbrt.org/fileformat-v3.html).

The Pbrt file format is an extensible file format for a plugin-based system.
Representing the format directly allows for best fidelity but pushed the burden
of interpreting standard plugins to the use. Yocto/ModelIO takes a different
approach and translates camera, shapes, materials, textures and lights from
Pbrt plugins to a common representation that presents users a simpler and
more uniform scene representation.

Yocto/ModelIO represents Pbrt data with the `pbrt_scene` struct.
Pbrt models are defined as collections of cameras, instanced shapes, materials,
texture and environments. Pbrt cameras are translate into a thin-len
approximations. Pbrt materials are translated to a material representation
similar to the Disney BSDF. Pbrt textures are either interpreted ion the fly or
defined by a image file. Pbrt area lights are translated to either emissive
materials and environments.

The Pbrt model is defined as an array of objects of the types defined above.
Pbrt objects are pointers owned by the main `pbrt_scene`.
Objects properties can be read and written directly from the model data,
and are documented in the header file for now.
Yocto/ModelIO does not currently provide functions to read and write Pbrt
shapes with a simpler interface than accessing data directly.

In general, Pbrt support is still experimental even if the library can
parse most Pbrt files. The objects properties documentations are for now
stored in the header file.

```cpp
auto pbrt = new pbrt_scene{...};            // obj model buffer
for(auto shape : pbrt->shapes)              // access shapes
  print_info(shape->name);                  // access shape properties
for(auto material : pbrt->material)         // access materials
  print_info(material->diffuse);            // access material properties
for(auto material : pbrt->material)         // access materials
  print_info(material->color_tex);          // access material textures
for(auto camera : pbrt->cameras)            // access cameras [extension]
  print_info(camera->frame);                // access camera properties
for(auto environment : pbrt->environments)  // access environments [extension]
  print_info(environment->emission);        // access environment properties
```

Use `load_pbrt(filename, pbrt, error)` to load Pbrt files and
`save_pbrt(filename, pbrt, error)` to save them.  
Both loading and saving take a filename, a pointer to a Pbrt model,
and returns whether or not the file was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.

```cpp
auto pbrt = new pbrt_scene{};           // obj model buffer
auto error = string{};                  // error buffer
if(!load_pbrt(filename, pbrt, error))   // load obj
  print_error(error);                   // check and print error
if(!save_pbrt(filename, pbrt, error))   // save obj
  print_error(error);                   // check and print error
```
