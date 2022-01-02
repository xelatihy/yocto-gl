# Yocto/ModelIO: Serialization for Obj, Ply and Stl models

Yocto/ModelIO is a collection of utilities for loading and saving scenes
and meshes in Ply, Obj, and Stl formats.
Yocto/ModelIO is implemented in `yocto_modelio.h` and `yocto_modelio.cpp`,
and depends on `fast_float.h`.

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

Use `ok = load_ply(filename, ply, error)` to load Ply files 
and `ok = save_ply(filename, ply, error)` to save them.  
Both loading and saving take a filename, a reference to a model, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.

```cpp
auto ply = ply_model{};              // ply model
auto error = string{};               // error string
if (!load_ply(filename, ply, error)) // load ply
  handle_error(error); 
if (!save_ply(filename, ply, error)) // save ply
  handle_error(error); 
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
auto ply = ply_model{};                 // ply model buffer
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
`get_triangles(ply, triangles)`, and `get_quads(ply, quads)`, to read
points, lines, triangles and quads. Note that since Ply support arbitrary
polygons and polylines, these functions tesselate the Ply polygons into
the desired element type, for now using a simple fan-like algorithm that
works only for convex elements. Use `has_quads(ply)` to check whether
the Ply data has quads, or use `get_faces(ply, triangles, quads)` to
triangles or quads.

```cpp
auto ply = ply_model{};                 // ply model buffer
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
auto ply = ply_model{};               // ply model buffer

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

auto error = string{};            // error buffer
save_ply(filename, ply, error);   // save ply
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
auto ply = ply_model{};                 // ply model buffer

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

Yocto/ModelIO represents Obj data with the `obj_model` struct.
Obj models are defined as collections of shapes and materials.
Obj shapes use a face-varying representation that has vertex positions,
normals and texture coordinates, with their their own topology.
In a shape, each face is tagged with a material used for that face.
Yocto/ModelIO provides direct access to these tagged shapes by inspecting the
`obj_shape` properties.

In Yocto/Obj, materials are represented by the `obj_material` struct.
Each material is a collection of values and textures that specify the material
lobes, like emission, diffuse, specular, reflection, etc. Each value
has a corresponding texture stored as an index to the texture array.

Yocto/ModelIO defines two extensions to the Obj file format.
Yocto/ModelIO adds another file, namely `.obx`, that stores cameras and
environment maps, respectively as `obj_camera` and `obj_environment`. T
his addition makes the extended Obj format capable of storing full scenes.

The Obj model is defined as an array of objects of the types defined above.
Obj objects are pointers owned by the main `obj_model`.
Objects properties can be read and written directly from the model data,
and are documented in the header file for now.
For shapes, Yocto/ModelIO provides several functions to read and write Obj
shapes, with a simpler interface than accessing data directly.

```cpp
auto obj = obj_model{...};                // obj model buffer
for(auto shape : obj.shapes)              // access shapes
  print_info(shape.name);                 // access shape properties
for(auto material : obj.material)         // access materials
  print_info(material.diffuse);           // access material properties
for(auto material : obj.material)         // access materials
  print_info(material.diffuse_tex);       // access material textures
for(auto camera : obj.cameras)            // access cameras [extension]
  print_info(camera.frame);               // access camera properties
for(auto environment : obj.environments)  // access environments [extension]
  print_info(environment.emission);       // access environment properties
```

Use `ok = load_obj(filename, obj, error)` to load Obj files
and `ok = save_obj(filename, obj, error)` to save them.
Both loading and saving take a filename, a reference to a model, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.

Obj is a face-varying file format, while most applications handle only
indexed meshes. The loading function takes an optional that specify whether
to load as face-varying or convert to indexed meshes, which is the default.

```cpp
auto obj = obj_model{};              // obj model
auto error = string{};               // error string
if (!load_obj(filename, obj, error)) // load obj
  handle_error(error); 
if (!save_obj(filename, obj, error)) // save obj
  handle_error(error); 
```

It is common in graphics to use Obj file to store single meshes. Yocto/Obj
supports this modality by providing specialized loading and saving functions
that take references to shapes as parameters.
Use `ok = load_obj(filename, shape, error)` to load Obj 
shapes and `ok = save_obj(filename, shape, error)` to save them. The loading function takes
an optional that specify whether to load as face-varying or convert to
indexed meshes, which is the default.

```cpp
auto oshape = obj_shape{};              // obj shape
auto error = string{};                  // error string
if (!load_obj(filename, oshape, error)) // load obj shape
  handle_error(error); 
if (!save_obj(filename, oshape, error)) // save obj shape
  handle_error(error); 
```

## Obj reading

Obj is a face-varying format and that geometry representation is maintained
in `obj_shape`. Yocto/ModelIO provides easier accessed to Obj shape data,
both as indexed meshes and as face-varying meshes.

Use `get_positions(shape, positions)`, `get_normals(shape, normals)`,
`get_texcoords(shape, texcoords, flipv)`, to read positions,
normals, texcoords, if present. Texture coordinates can be
optionally flipped vertically. For shape elements, use
`get_points(shape, points, materials)`, `get_lines(shape, lines, materials)`,
`get_triangles(shape, triangles, materials)`, and
`get_quads(shape, quads, materials)`, to read points, lines, triangles and quads.
Note that since Obj support arbitrary polygons and polylines, these functions
tesselate the Obj polygons into the desired element type, for now using a
simple fan-like algorithm that works only for convex elements.
Use `has_quads(shape)` to check whether the Obj data has quads,
or use `get_faces(shape, triangles, quads, materials)` to triangles or quads.

In some cases, it may be desireable to extract the shape elements corresponding
to a single material, for example for use in renderers that support a single
shader per shape. To filtering elements by material index use
`get_points(shape, material, points)`, `get_lines(shape, material, lines)`,
`get_triangles(shape, material, triangles)`, and `get_quads(shape, material, quads)`.
that are overrides of the previous `get_<element>(...)` functions,
but differ in that they take a material is as input, instead of returning
materials tags.

```cpp
auto obj = obj_model{};                // obj model buffer
auto error = string{};                 // error buffer
load_obj(filename, obj, error);        // load obj
auto& shape = obj.shapes.front();      // get shape

auto positions = vector<vec3f>{};      // vertex properties
get_positions(shape, positions);
auto normals   = vector<vec3f>{};
get_normals(shape, normals);
auto texcoords = vector<vec2f>{};
get_texcoords(shape, texcoords);
auto triangles = vector<vec3i>{};      // element data
auto quads     = vector<vec4i>{};
auto materials = vector<int>{};        // per-face material ids
if(has_quads(shape)) {
  get_triangles(shape, triangles, materials); // read as triangles
} else {
  get_quads(shape, quads, materials);   // read as quads
}

auto material_id = 0;                     // material id to extract to
if(has_quads(shape)) {
  get_triangles(shape, material_id, triangles); // read as triangles
} else {
  get_quads(shape, material_id, quads);   // read as quads
}
```

Yocto/ModelIO supports also reading Obj shapes as face-varying quads
with `get_fvquads(...)`. In this case, the model should be loaded as
face-varying.

```cpp
auto obj =  obj_model{};                // obj model buffer
auto error = string{};                  // error buffer
load_obj(filename, obj, error, true);   // load obj as face-varying
auto& shape = obj.shapes.front();       // get shape

auto positions = vector<vec3f>{};       // vertex properties
get_positions(shape, positions);
auto normals   = vector<vec3f>{};
get_normals(shape, normals);
auto texcoords = vector<vec2f>{};
get_texcoords(shape, texcoords);
auto quadspos  = vector<vec4i>{};       // face-varying element data
auto quadsnorm = vector<vec4i>{};
auto quadsuv   = vector<vec4i>{};
auto materials = vector<int>{};         // per-face material ids
get_fvquads(shape,                      // read as face-varying quads
  quadspos, quadsnorm, quadsuv, materias);
```

## Obj writing

To save an Obj, create a scene and add objects to it by manipulating the
objects' arrays directly. For all objects, set the objects' properties directly.

For shapes, Yocto/ModelIO defines convenience functions that take either
indexed mesh or face-varying meshes as input and create the appropriate
Obj shape elements. Use `add_positions(shape, positions)`,
`add_normals(shape, normals)` and `add_texcoords(shape, texcoords, flipv)`
to add vertex properties to a shape.
Use `add_triangles(shape, triangles, material, has_normals, has_texcoords)`,
`add_quads(shape, quads, material, has_normals, has_texcoords)`,
`add_lines(shape, lines, material, has_normals, has_texcoords)`,
`add_points(shape, points, material, has_normals, has_texcoords)`
to add triangles, quads, lines or points to shapes respectively.
Material data is only represented as tags and can be left
empty if only one material is used. To set material names
use `set_materials(shape, materials)`.
To add face-varying shapes, use `add_fvquads(...)`.

```cpp
auto obj =  obj_model{};                 // obj model

auto camera = obj_camera{};              // init camera
camera.name = "camera";                  // set camera name
camera.frame = identity3x4f;             // set camera properties
obj.cameras.push_back(camera);           // add camera
auto environment = add_environment(obj); // init environment
environment.name = "environment";        // set environment name
environment.emission = {1,1,1};          // set environment properties
obj.environments.push_back(environment); // add environment
auto material = add_material(obj);       // init material
material.name = "material";              // set material name
material.diffuse = {1,0,0};              // set material properties
obj.materials.push_back(material);       // add material

auto triangles = vector<vec3i>{...};     // element data
auto positions = vector<vec3f>{...};     // vertex properties
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>{...};

auto shape = obj_shape{};                // init shape
shape.name = "shape";                    // set shape name
add_triangles(shape, triangles, 0);      // set shape geometry
add_positions(shape, positions);
add_normals(shape, normals);
add_texcoords(shape, texcoords);

auto error = string{};                  // error buffer
save_obj(filename, obj, error);         // save obj
```

And similarly for face-varying shapes.

```cpp
auto obj = obj_model{};                 // obj model buffer

auto quadspos  = vector<vec4i>{...};    // face-varying element data
auto quadsnorm = vector<vec4i>{...};
auto quadsuv   = vector<vec4i>{...};
auto positions = vector<vec3f>{...};    // vertex properties
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>....{};

auto shape = obj_shape{};               // add shape
shape.name = "shape";                   // set shape name
add_fvquads(shape, quadspos,            // set shape geometry
  quadsnorm, quadstexcoord, 0);
add_positions(shape, positions);
add_normals(shape, normals);
add_texcoords(shape, texcoords);

auto error = string{};                  // error buffer
save_obj(filename, obj, error);         // save obj
```

## Stl models

The Stl file format is a file format used to serialize meshes.
To use this library is helpful to understand the basic of the Stl
file format for example from the
[Stl Wikipedia page](<https://en.wikipedia.org/wiki/STL_(file_format)>).

Yocto/ModelIO represents Stl data with the `stl_model` struct.
Stl models are defined as collections of shapes, stored in `stl_shape`.
Shapes are defined by an array of vertex positions, an array of triangle
indices, and an array of triangle normals.

Use `ok = load_stl(filename, stl, error)` to load Stl files
and `ok = save_stl(filename, stl, error)` to save them.

```cpp
auto stl = stl_model{};              // stl model
auto error = string{};               // error string
if (!load_stl(filename, stl, error)) // load stl
  handle_error(error); 
if (!save_stl(filename, stl, error)) // save stl
  handle_error(error); 
```

## Stl reading

You can access Stl data directly from the shapes, without any further complexity.

```cpp
auto stl = stl_model{};                // stl model buffer
load_obj(filename, stl);               // load stl
auto& shape = stl.shapes.front();      // get shape

auto positions = shape.positions;      // vertex properties
auto triangles = shape.triangles;      // element data
```

## Stl writing

Similarly, you can create an Stl model by acting directly on the data.

```cpp
auto stl = stl_model{};                // stl model buffer
auto shape = stl_shape{};              // create shape
shape.positions = vector<vec3f>{...};
shape.triangles = vector<vec3i>{...};
stl.shapes.push_back(shape);           // add shape
load_obj(filename, stl);               // save stl
```
