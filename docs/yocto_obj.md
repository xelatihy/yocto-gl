# Yocto/Obj

Wavefront OBJ/MTL loader and writer with support for points,
lines, triangles and general polygons and all materials properties.
Contains also a few extensions to easily create demos such as per-vertex
color and radius, cameras, environment maps and instances.
Can use either a low-level OBJ representation or a high level flattened
representation.

Both in reading and writing, OBJ has no clear convention on the orientation
of textures Y axis. So in many cases textures appears flipped. To handle
that, use the option to flip textures coordinates on either saving or
loading. By default texture coordinates are flipped since this seems
the convention found on test cases collected on the web. The value Tr
has similar problems, since its relation to opacity is software specific.
Again we let the user chose the convension and set the default to the
one found on the web.

In the high level interface, shapes are indexed meshes and are described
by arrays of vertex indices for points/lines/triangles and arrays for vertex
positions, normals, texcoords, color and radius. The latter two as
extensions. Since OBJ is a complex formats that does not match well with
current GPU rendering / path tracing algorithms, we adopt a simplification
similar to other single file libraries:
1. vertex indices are unique, as in OpenGL and al standard indexed triangle
  meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
  vertex dusplication happens thought for same triplets
2. we split shapes on changes to groups and materials, instead of keeping
  per-face group/material data; this makes the data usable right away in
  a GPU viewer; this is not a major limitation if we accept the previous
  point that already changes shapes topology.

This library depends in yocto_math.h. Texture loading depends on
yocto_image. If the texture loading dependency is not desired, it can be
disabled by defining YOBJ_NO_IMAGE before including this file.


## Usage Of High-Level Interface

1. load a scene with `load_scene()`
2. look at the `scene` data structures for access to individual elements
3. can also manipulate the scene by adding missing data with `add_XXX()`
   functions
4. since OBJ does natively support mesh transfotms, which we support with
   instances, use `flatten_instaces()` or `add_instances()` to go back
   and fourth
5. use `save_scene()` ti write the data to disk

## Usage Of Low-Level Interface

1. load a obj data with `load_obj()`; can load also textues
2. look at the `obj_XXX` data structures for access to individual elements
3. use obj back to disk with `save_obj()`; can also save textures
4. conversion from low- to -high-level data structures with
   `scene_to_obj()` and `obj_to_scene()`


## History

- v 0.30: support for smoothing groups
- v 0.29: use reference interface for textures
- v 0.28: add function to split meshes into single shapes
- v 0.27: explicit transforms
- v 0.26: added interpreting of illum in scene conversions
- v 0.25: added convention for Tr
- v 0.24: remove exception from code and add explicit error handling
- v 0.23: texture have always 4 channels
- v 0.22: change variable names for compilation on gcc
- v 0.21: bug fixes
- v 0.20: use yocto_math in the interface and remove inline compilation
- v 0.19: add missing bounding box computation and missing data functions
- v 0.18: prioritize high-level interface
- v 0.17: name cleanup in both interface to better align with glTF
- v 0.16: change flattened data structure to use pointers
- v 0.15: added unknown properties std::map for materials
- v 0.14: added extension to store tetrahedral meshes
- v 0.13: started adding physics extension to materials
- v 0.12: change texture loading by flipping uvs rather than images
- v 0.11: use yocto_image for texture handling.
- v 0.10: switch to .h/.cpp pair
- v 0.9: bug fixes and optionally texture skipping
- v 0.8: high level interface uses grouping
- v 0.7: doxygen comments
- v 0.6: bug fixes
- v 0.5: removed options to force image formats (image library not reliable)
- v 0.4: [major API change] move to modern C++ interface
- v 0.3: new API internals and C++ interface
- v 0.2: removal of C interface
- v 0.1: C++ implementation
- v 0.0: initial release in C99

## Namespace yobj

Reading and Writing support for Wavefront OBJ.

### Typedef template  <typename T \>
property_map

~~~ .cpp
template <typename T>
using property_map = std::map<std::string, std::vector<T>>;
~~~

Property map

### Struct texture

~~~ .cpp
struct texture {
    std::string path;
    ym::image4b ldr;
    ym::image4f hdr;
    int width() const; 
    int height() const; 
}
~~~

Scene Texture

- Members:
    - path:      path
    - ldr:      if loaded, ldr image
    - hdr:      if loaded, hdr image
    - width():      get texture width
    - height():      get texture height


### Struct texture_info

~~~ .cpp
struct texture_info {
    bool clamp = false;
    float bump_scale = 1;
    property_map<std::string> unknown_props;
}
~~~

Scene Texture Additional Information

- Members:
    - clamp:      clamping
    - bump_scale:      bump scale
    - unknown_props:      unknown string props


### Struct material

~~~ .cpp
struct material {
    std::string name;
    ym::vec3f ke = {0, 0, 0};
    ym::vec3f kd = {0, 0, 0};
    ym::vec3f ks = {0, 0, 0};
    ym::vec3f kt = {0, 0, 0};
    float rs = 0.0001;
    float opacity = 1;
    texture* ke_txt = nullptr;
    texture* kd_txt = nullptr;
    texture* ks_txt = nullptr;
    texture* kt_txt = nullptr;
    texture* rs_txt = nullptr;
    texture* op_txt = nullptr;
    texture* bump_txt = nullptr;
    texture* disp_txt = nullptr;
    texture* norm_txt = nullptr;
    texture_info ke_txt_info = {};
    texture_info kd_txt_info = {};
    texture_info ks_txt_info = {};
    texture_info kt_txt_info = {};
    texture_info rs_txt_info = {};
    texture_info bump_txt_info = {};
    texture_info disp_txt_info = {};
    texture_info norm_txt_info = {};
    property_map<std::string> unknown_props;
}
~~~

Scene Material

- Members:
    - name:      material name
    - ke:      emission color
    - kd:      diffuse color
    - ks:      specular color
    - kt:      transmission color
    - rs:      roughness
    - opacity:      opacity
    - ke_txt:      emission texture
    - kd_txt:      diffuse texture
    - ks_txt:      specular texture
    - kt_txt:      transmission texture
    - rs_txt:      roughness texture
    - op_txt:      opacity texture
    - bump_txt:      bump map texture (heighfield)
    - disp_txt:      displacement map texture (heighfield)
    - norm_txt:      normal texture
    - ke_txt_info:      emission texture
    - kd_txt_info:      diffuse texture
    - ks_txt_info:      specular texture
    - kt_txt_info:      transmission texture
    - rs_txt_info:      roughness texture
    - bump_txt_info:      bump map texture (heighfield)
    - disp_txt_info:      displacement map texture (heighfield)
    - norm_txt_info:      normal texture
    - unknown_props:      unknown string props


### Struct shape

~~~ .cpp
struct shape {
    std::string name = "";
    material* mat = nullptr;
    std::vector<int> points;
    std::vector<ym::vec2i> lines;
    std::vector<ym::vec3i> triangles;
    std::vector<ym::vec4i> tetras;
    std::vector<ym::vec3f> pos;
    std::vector<ym::vec3f> norm;
    std::vector<ym::vec2f> texcoord;
    std::vector<ym::vec4f> color;
    std::vector<float> radius;
    std::vector<ym::vec4f> tangsp;
}
~~~

Shape. May contain only one of the points/lines/triangles.

- Members:
    - name:      name of the group that enclosed it
    - mat:      material
    - points:      points
    - lines:      lines
    - triangles:      triangles
    - tetras:      tetrahedrons
    - pos:      per-vertex position (3 float)
    - norm:      per-vertex normals (3 float)
    - texcoord:      per-vertex texcoord (2 float)
    - color:      [extension] per-vertex color (4 float)
    - radius:      [extension] per-vertex radius (1 float)
    - tangsp:      [extension] per-vertex tangent space (4 float)


### Struct mesh

~~~ .cpp
struct mesh {
    std::vector<shape*> shapes;
    ~mesh(); 
}
~~~

Mesh

- Members:
    - shapes:      primitives
    - ~mesh():      cleanup


### Struct instance

~~~ .cpp
struct instance {
    ym::vec3f translation = {0, 0, 0};
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::vec3f scale = {1, 1, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    mesh* msh = nullptr;
    ym::mat4f xform() const; 
}
~~~

Mesh instance.

- Members:
    - translation:      translation
    - rotation:      rotation
    - scale:      scale
    - matrix:      generic transform matrix
    - msh:      mesh instances
    - xform():      xform


### Struct camera

~~~ .cpp
struct camera {
    std::string name;
    ym::vec3f translation = {0, 0, 0};
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    bool ortho = false;
    float yfov = 2;
    float aspect = 16.0f / 9.0f;
    float focus = 1;
    float aperture = 0;
    ym::mat4f xform() const; 
}
~~~

Scene Camera

- Members:
    - name:      name
    - translation:      translation
    - rotation:      rotation
    - matrix:      generic transform matrix
    - ortho:      ortho cam
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - focus:      focus distance
    - aperture:      lens aperture
    - xform():      xform


### Struct environment

~~~ .cpp
struct environment {
    std::string name;
    material* mat = nullptr;
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    ym::mat4f xform() const; 
}
~~~

Envinonment map

- Members:
    - name:      name
    - mat:      index of material in material array
    - rotation:      rotation
    - matrix:      generic transform matrix
    - xform():      xform


### Struct scene

~~~ .cpp
struct scene {
    std::vector<mesh*> meshes;
    std::vector<instance*> instances;
    std::vector<material*> materials;
    std::vector<texture*> textures;
    std::vector<camera*> cameras;
    std::vector<environment*> environments;
    ~scene(); 
}
~~~

Scene

- Members:
    - meshes:      shape array
    - instances:      instance array
    - materials:      material array
    - textures:      texture array
    - cameras:      camera array
    - environments:      environment array
    - ~scene():      cleanup


### Function load_scene()

~~~ .cpp
scene* load_scene(const std::string& filename, bool load_textures,
    bool flip_texcoord = true, bool skip_missing = true, bool flip_tr = true,
    std::string* err = nullptr);
~~~

Load scene

- Parameters:
    - filename: filename
    - load_textures: whether to load textures (default to false)
    - flip_texcoord: whether to flip the v coordinate
    - skip_missing: skip missing textures
    - flip_tr: whether to flip tr
    - err: if set, store error message on error
- Returns:
    - scene (nullptr on error)

### Function save_scene()

~~~ .cpp
bool save_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool flip_texcoord = true, bool flip_tr = true,
    std::string* err = nullptr);
~~~

Save scene

- Parameters:
    - filename: filename
    - scn: scene data to save
    - save_textures: whether to save textures (default to false)
    - flip_texcoord: whether to flip the v coordinate
    - flip_tr: whether to flip tr
    - err: if set, store error message on error
- Returns:
    - whether an error occurred

### Function load_textures()

~~~ .cpp
bool load_textures(scene* scn, const std::string& dirname,
    bool skip_missing = false, std::string* err = nullptr);
~~~

Loads textures for an scene.

- Parameters:
    - scn: scene to load textures into
    - dirname: base directory name for texture files
    - skip_missing: whether to skip missing textures or stops with error
    - err: if set, store error message on error
- Returns:
    - whether an error occurred

### Function save_textures()

~~~ .cpp
bool save_textures(const scene* scn, const std::string& dirname,
    bool skip_missing = false, std::string* err = nullptr);
~~~

Saves textures for an scene.

- Parameters:
    - scn: scene to write textures from
    - dirname: base directory name for texture files
    - skip_missing: whether to skip missing textures or stops with error
    - err: if set, store error message on error
- Returns:
    - whether an error occurred

### Function compute_scene_bounds()

~~~ .cpp
ym::bbox3f compute_scene_bounds(const scene* scn);
~~~

Computes a scene bounding box

### Function add_normals()

~~~ .cpp
void add_normals(scene* scn);
~~~

Add missing data to the scene.

### Function add_radius()

~~~ .cpp
void add_radius(scene* scn, float radius);
~~~

Add missing data to the scene.

### Function add_tangent_space()

~~~ .cpp
void add_tangent_space(scene* scn);
~~~

Add missing data to the scene.

### Function add_texture_data()

~~~ .cpp
void add_texture_data(scene* scn);
~~~

Add missing data to the scene.

### Function add_instances()

~~~ .cpp
void add_instances(scene* scn);
~~~

Add missing data to the scene.

### Function add_names()

~~~ .cpp
void add_names(scene* scn);
~~~

Add missing data to the scene.

### Function add_default_camera()

~~~ .cpp
void add_default_camera(scene* scn);
~~~

Add a default camera that views the entire scene.

### Function flatten_instances()

~~~ .cpp
void flatten_instances(scene* scn);
~~~

Flatten scene instances into separate meshes.

### Function split_shapes()

~~~ .cpp
void split_shapes(scene* scn);
~~~

Split meshes into single shapes

### Struct obj_vertex

~~~ .cpp
struct obj_vertex {
    int pos;
    int texcoord;
    int norm;
    int color;
    int radius;
    obj_vertex(int pos = -1, int texcoord = -1, int norm = -1, int color = -1, int radius = -1); 
}
~~~

Face vertex

- Members:
    - pos:      position
    - texcoord:      texcoord
    - norm:      normal
    - color:      color [extension]
    - radius:      radius [extension]
    - obj_vertex():      Constructor (copies members initializing missing ones to -1)


### Enum obj_element_type : uint16_t

~~~ .cpp
enum struct obj_element_type : uint16_t {
    point = 1,
    line = 2,
    face = 3,
    tetra = 4,
}
~~~

element type

- Values:
    - point:      lists of points
    - line:      polylines
    - face:      polygon faces
    - tetra:      tetrahedrons


### Struct obj_element

~~~ .cpp
struct obj_element {
    uint32_t start;
    obj_element_type type;
    uint16_t size;
}
~~~

Element vertex indices

- Members:
    - start:      starting vertex index
    - type:      element type
    - size:      number of vertices


### Struct obj_group

~~~ .cpp
struct obj_group {
    std::string matname;
    std::string groupname;
    bool smoothing = true;
    std::vector<obj_vertex> verts;
    std::vector<obj_element> elems;
}
~~~

Element group

- Members:
    - matname:      material name
    - groupname:      group name
    - smoothing:      smoothing
    - verts:      element vertices
    - elems:      element faces


### Struct obj_object

~~~ .cpp
struct obj_object {
    std::string name;
    std::vector<obj_group> groups;
}
~~~

Obj object

- Members:
    - name:      object name
    - groups:      element groups


### Struct obj_material

~~~ .cpp
struct obj_material {
    std::string name;
    int illum = 0;
    ym::vec3f ke = {0, 0, 0};
    ym::vec3f ka = {0, 0, 0};
    ym::vec3f kd = {0, 0, 0};
    ym::vec3f ks = {0, 0, 0};
    ym::vec3f kr = {0, 0, 0};
    ym::vec3f kt = {0, 0, 0};
    float ns = 1;
    float ior = 1;
    float op = 1;
    std::string ke_txt;
    std::string ka_txt;
    std::string kd_txt;
    std::string ks_txt;
    std::string kr_txt;
    std::string kt_txt;
    std::string ns_txt;
    std::string op_txt;
    std::string ior_txt;
    std::string bump_txt;
    std::string disp_txt;
    std::string norm_txt;
    property_map<std::string> ke_txt_info = {};
    property_map<std::string> ka_txt_info = {};
    property_map<std::string> kd_txt_info = {};
    property_map<std::string> ks_txt_info = {};
    property_map<std::string> kr_txt_info = {};
    property_map<std::string> kt_txt_info = {};
    property_map<std::string> ns_txt_info = {};
    property_map<std::string> op_txt_info = {};
    property_map<std::string> ior_txt_info = {};
    property_map<std::string> bump_txt_info = {};
    property_map<std::string> disp_txt_info = {};
    property_map<std::string> norm_txt_info = {};
    property_map<std::string> unknown_props;
}
~~~

OBJ material

- Members:
    - name:      material name
    - illum:      MTL illum mode
    - ke:      emission color
    - ka:      ambient color
    - kd:      diffuse color
    - ks:      specular color
    - kr:      reflection color
    - kt:      transmision color
    - ns:      phong exponent for ks
    - ior:      index of refraction
    - op:      opacity
    - ke_txt:      emission texture
    - ka_txt:      ambient texture
    - kd_txt:      diffuse texture
    - ks_txt:      specular texture
    - kr_txt:      reflection texture
    - kt_txt:      transmission texture
    - ns_txt:      specular exponent texture
    - op_txt:      opacity texture
    - ior_txt:      index of refraction
    - bump_txt:      bump map texture (heighfield)
    - disp_txt:      displacement map texture (heighfield)
    - norm_txt:      normal map texture
    - ke_txt_info:      emission texture
    - ka_txt_info:      ambient texture
    - kd_txt_info:      diffuse texture
    - ks_txt_info:      specular texture
    - kr_txt_info:      reflection texture
    - kt_txt_info:      transmission texture
    - ns_txt_info:      specular exponent texture
    - op_txt_info:      opacity texture
    - ior_txt_info:      index of refraction
    - bump_txt_info:      bump map texture (heighfield)
    - disp_txt_info:      displacement map texture (heighfield)
    - norm_txt_info:      normal texture
    - unknown_props:      unknown string props


### Struct obj_camera

~~~ .cpp
struct obj_camera {
    std::string name;
    ym::vec3f translation = {0, 0, 0};
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::vec3f scale = {1, 1, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    bool ortho = false;
    float yfov = 2 * std::atan(0.5f);
    float aspect = 16.0f / 9.0f;
    float aperture = 0;
    float focus = 1;
}
~~~

Camera [extension]

- Members:
    - name:      camera name
    - translation:      translation
    - rotation:      rotation
    - scale:      scale
    - matrix:      generic transform matrix
    - ortho:      orthografic camera
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - aperture:      lens aperture
    - focus:      focus distance


### Struct obj_environment

~~~ .cpp
struct obj_environment {
    std::string name;
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    std::string matname;
}
~~~

Environment [extension]

- Members:
    - name:      environment name
    - rotation:      rotation
    - matrix:      generic transform matrix
    - matname:      material name


### Struct obj_instance

~~~ .cpp
struct obj_instance {
    std::string name;
    ym::vec3f translation = {0, 0, 0};
    ym::quat4f rotation = {0, 0, 0, 1};
    ym::vec3f scale = {1, 1, 1};
    ym::mat4f matrix = ym::identity_mat4f;
    std::string meshname;
}
~~~

Instance [extension]

- Members:
    - name:      instance name
    - translation:      translation
    - rotation:      rotation
    - scale:      scale
    - matrix:      generic transform matrix
    - meshname:      object name


### Struct obj

~~~ .cpp
struct obj {
    std::vector<ym::vec3f> pos;
    std::vector<ym::vec3f> norm;
    std::vector<ym::vec2f> texcoord;
    std::vector<ym::vec4f> color;
    std::vector<float> radius;
    std::vector<obj_object> objects;
    std::vector<obj_material> materials;
    std::vector<obj_camera> cameras;
    std::vector<obj_environment> environments;
    std::vector<obj_instance> instances;
}
~~~

OBJ asset

- Members:
    - pos:      vertex positions
    - norm:      vertex normals
    - texcoord:      vertex texcoord
    - color:      vertex color [extension]
    - radius:      vertex radius [extension]
    - objects:      objects
    - materials:      materials
    - cameras:      cameras [extension]
    - environments:      env maps [extension]
    - instances:      instances [extension]


### Function load_obj()

~~~ .cpp
obj* load_obj(const std::string& filename, bool flip_texcoord = true,
    bool flip_tr = true, std::string* err = nullptr);
~~~

Load OBJ

- Parameters:
    - filename: filename
    - flip_texcoord: whether to flip the v coordinate
    - flip_tr: whether to flip the Tr value
    - err: if set, store error message on error
- Return:
    - obj (nullptr on error)

### Function load_mtl()

~~~ .cpp
std::vector<obj_material> load_mtl(const std::string& filename,
    bool flip_tr = true, std::string* err = nullptr);
~~~

Load MTL

- Parameters:
    - filename: filename
    - flip_tr: whether to flip the Tr value
    - err: if set, store error message on error
- Return:
    - loaded materials (empty on error)

### Function save_obj()

~~~ .cpp
bool save_obj(const std::string& filename, const obj* asset,
    bool flip_texcoord = true, bool flip_tr = true, std::string* err = nullptr);
~~~

Save OBJ

- Parameters:
    - filename: filename
    - asset: obj data to save
    - flip_texcoord: whether to flip the v coordinate
    - flip_tr: whether to flip the Tr value
    - err: if set, store error message on error
- Returns:
    - whether an error occurred

### Function save_mtl()

~~~ .cpp
bool save_mtl(const std::string& filename,
    const std::vector<obj_material>& materials, bool flip_tr = true,
    std::string* err = nullptr);
~~~

Save MTL (@deprecated interface)

### Function obj_to_scene()

~~~ .cpp
scene* obj_to_scene(const obj* obj);
~~~

Converts an OBJ into a scene.

- Parameters:
    - obj: obj to be flattened

### Function scene_to_obj()

~~~ .cpp
obj* scene_to_obj(const scene* scn);
~~~

Save a scene in an OBJ file.

- Parameters:
    - scn: scene to convert

