# Yocto/ModelIO: Serialization for Obj, Ply and Pbrt models

Yocto/ModelIO is a collection of utilities for loading and saving scenes
and meshes in Ply, Obj and Pbrt formats.
Yocto/ModelIO is implemented in `yocto_modelio.h` and `yocto_modelio.cpp`.

## Loading and saving Ply files

The Ply file format is a generic file format used to serialize meshes and point
clouds. To use this library is helpful to understand the basic of the Ply
file format for example from [this Wikipedia page](<https://en.wikipedia.org/wiki/PLY_(file_format)>).

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

Yocto/Shape defines convenience functions to add the most used properties
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

## Obj IO

// Load and save obj
bool load_obj(const string& filename, obj_model* obj, string& error,
bool geom_only = false, bool split_elements = true,
bool split_materials = false);
bool save_obj(const string& filename, obj_model* obj, string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_triangles(const obj_shape* shape, vector<vec3i>& triangles,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<obj_material*>& materials, vector<int>& ematerials,
bool flip_texcoord = false);
void get_quads(const obj_shape* shape, vector<vec4i>& quads,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<obj_material*>& materials, vector<int>& ematerials,
bool flip_texcoord = false);
void get_lines(const obj_shape* shape, vector<vec2i>& lines,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<obj_material*>& materials, vector<int>& ematerials,
bool flip_texcoord = false);
void get_points(const obj_shape* shape, vector<int>& points,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<obj_material*>& materials, vector<int>& ematerials,
bool flip_texcoord = false);
void get_fvquads(const obj_shape* shape, vector<vec4i>& quadspos,
vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<obj_material*>& materials, vector<int>& ematerials,
bool flip_texcoord = false);
bool has_quads(obj_shape\* shape);

// Get obj shape by extracting the elements beloing to only one material.
void get_triangles(const obj_shape* shape, int material,
vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_quads(const obj_shape* shape, int material, vector<vec4i>& quads,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
bool flip_texcoord = false);
void get_lines(const obj_shape* shape, int material, vector<vec2i>& lines,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
bool flip_texcoord = false);
void get_points(const obj_shape* shape, int material, vector<int>& points,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
bool flip_texcoord = false);

// Create OBJ
obj_camera* add_camera(obj_model* obj);
obj_material* add_material(obj_model* obj);
obj_environment* add_environment(obj_model* obj);
obj_shape* add_shape(obj_model* obj);

// Add obj shape
void set_triangles(obj_shape* shape, const vector<vec3i>& triangles,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
bool flip_texcoord = false);
void set_quads(obj_shape* shape, const vector<vec4i>& quads,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
bool flip_texcoord = false);
void set_lines(obj_shape* shape, const vector<vec2i>& lines,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
bool flip_texcoord = false);
void set_points(obj_shape* shape, const vector<int>& points,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
bool flip_texcoord = false);
void set_fvquads(obj_shape* shape, const vector<vec4i>& quadspos,
const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
bool flip_texcoord = false);
void set_materials(obj_shape* shape, const vector<obj_material*>& materials);
void set_instances(obj_shape* shape, const vector<frame3f>& instances);

## Pbrt IO

// Load/save pbrt
bool load_pbrt(const string& filename, pbrt_model* pbrt, string& error);
bool save_pbrt(const string& filename, pbrt_model* pbrt, string& error,
bool ply_meshes = false);

// Create pbrt
pbrt_camera* add_camera(pbrt_model* pbrt);
pbrt_shape* add_shape(pbrt_model* pbrt);
pbrt_material* add_material(pbrt_model* pbrt);
pbrt_environment* add_environment(pbrt_model* pbrt);
pbrt_light* add_light(pbrt_model* pbrt);

```

```
