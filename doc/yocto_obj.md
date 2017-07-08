# Yocto/Obj

Wavefront OBJ/MTL loader and writer with support for points,
lines, triangles and general polygons and all materials properties.
Contains also a few extension to eqasily create demos such as per-vertex
color and radius, cameras and envmaps. Can use either a low-level OBJ
representation or a high level flattened representation.

Both in reading and writing, OBJ has no clear convention on the orientation
of textures Y axis. So in many cases textures appears flipped. To handle
that, use the option to flip textures coordinates on either saving or
loading. By default texture coordinates are flipped since this seems
the convention found on test cases collected on the web.

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


## Usage for reading

1. load an obj with load_obj()
    - loads an obj from disk including its associate mtl files
    - returns a parsed scene data structure described below
2. [LOW-LEVEL INTERFACE] access the data directly from the returned object
    - the data is documented below and matches the OBJ file structure
    exactly
3. [HIGH-LEVEL INTERFACE] optionally flatten the data as a more friendly
     representation where shapes are index meshes, supporting points, lines
     and triangle primitives, with flatten_obj()
    - the flattened data, documented below, can be use to draw directly on
      the GPU or in a raytracer
    - vertices are duplicated as needed to support GPU friendly access
    - optionally load textures data as float arrays with
      load_fl_textures()

## Usage for Writing

1. include this file (more compilation options below)
2. [LOW-LEVEL INTERFACE] fill an obj object with your scene data and save
   the obj/mtl pair with save_obj()
   ok = save_obj(filename, obj, error message)
3. [HIGH_LEVEL INTERFACE] create a flattened scene object and turn into an
   obj with unflatten_obj()


## History

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
    int width = 0;
    int height = 0;
    int ncomp = 0;
    std::vector<float> dataf;
    std::vector<unsigned char> datab;
}
~~~

Scene Texture

- Members:
    - path:      path
    - width:      if loaded, image width
    - height:      if loaded, image hieght
    - ncomp:      if loaded, number of component (1-4)
    - dataf:      if loaded, pixel data for HDRs
    - datab:      if loaded, pixel data for LDRs


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
    texture* bump_txt = nullptr;
    texture* disp_txt = nullptr;
    texture* norm_txt = nullptr;
    float stiffness = 0.0f;
    float density = 0.0f;
    property_map<std::string> str_props;
    property_map<int> int_props;
    property_map<float> flt_props;
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
    - bump_txt:      bump map texture (heighfield)
    - disp_txt:      displacement map texture (heighfield)
    - norm_txt:      normal texture
    - stiffness:      stiffness
    - density:      density
    - str_props:      unknown string props
    - int_props:      unknown int props
    - flt_props:      unknown float props


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
    ym::mat4f xform = ym::identity_mat4f;
    mesh* mesh = nullptr;
}
~~~

Mesh instance.

- Members:
    - xform:      transform
    - mesh:      primitives


### Struct camera

~~~ .cpp
struct camera {
    std::string name;
    ym::mat4f xform = ym::identity_mat4f;
    bool ortho = false;
    float yfov = 2;
    float aspect = 16.0f / 9.0f;
    float focus = 1;
    float aperture = 0;
}
~~~

Scene Camera

- Members:
    - name:      name
    - xform:      transform
    - ortho:      ortho cam
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - focus:      focus distance
    - aperture:      lens aperture


### Struct environment

~~~ .cpp
struct environment {
    std::string name;
    material* mat = nullptr;
    ym::mat4f xform = ym::identity_mat4f;
}
~~~

Envinonment map

- Members:
    - name:      name
    - mat:      index of material in material array
    - xform:      transform


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
scene* load_scene(
    const std::string& filename, bool load_textures, bool flip_texcoord = true);
~~~

Load scene

- Parameters:
    - filename: filename
    - load_textures: whether to load textures (default to false)
    - flip_texcoord: whether to flip the v coordinate
- Return:
    - scene

### Function save_scene()

~~~ .cpp
void save_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool flip_texcoord = true);
~~~

Save scene

- Parameters:
    - filename: filename
    - scn: scene data to save
    - save_textures: whether to save textures (default to false)
    - flip_texcoord: whether to flip the v coordinate

### Function load_textures()

~~~ .cpp
void load_textures(
    scene* scn, const std::string& dirname, bool skip_missing = false);
~~~

Loads textures for an scene.

- Parameters:
    - scn: scene to load textures into
    - dirname: base directory name for texture files
    - skip_missing: whether to skip missing textures or throw an expection

### Function save_textures()

~~~ .cpp
void save_textures(
    const scene* scn, const std::string& dirname, bool skip_missing = false);
~~~

Saves textures for an scene.

- Parameters:
    - scn: scene to write textures from
    - dirname: base directory name for texture files
    - skip_missing: whether to skip missing textures or throw an expection

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
    obj_element(uint32_t start_, etype type_, uint16_t size_); 
}
~~~

Element vertex indices

- Members:
    - start:      starting vertex index
    - type:      element type
    - size:      number of vertices
    - obj_element():      constructor


### Struct obj_group

~~~ .cpp
struct obj_group {
    std::string matname;
    std::string groupname;
    std::vector<obj_vertex> verts;
    std::vector<obj_element> elems;
}
~~~

Element group

- Members:
    - matname:      material name
    - groupname:      group name
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


### Struct mtl_material

~~~ .cpp
struct mtl_material {
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
    float stiffness = 0;
    float density = 0;
    property_map<std::string> str_props;
    property_map<int> int_props;
    property_map<float> flt_props;
}
~~~

MTL material

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
    - stiffness:      overall stiffness
    - density:      density
    - str_props:      unknown string props
    - int_props:      unknown int props
    - flt_props:      unknown float props


### Struct obj_camera

~~~ .cpp
struct obj_camera {
    std::string name;
    ym::mat4f xform = ym::identity_mat4f;
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
    - xform:      camera transform
    - ortho:      orthografic camera
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - aperture:      lens aperture
    - focus:      focus distance


### Struct obj_environment

~~~ .cpp
struct obj_environment {
    std::string name;
    ym::mat4f xform = ym::identity_mat4f;
    std::string matname;
}
~~~

Environment [extension]

- Members:
    - name:      environment name
    - xform:      transform
    - matname:      material name


### Struct obj_instance

~~~ .cpp
struct obj_instance {
    std::string name;
    ym::mat4f xform = ym::identity_mat4f;
    std::string meshname;
}
~~~

Instance [extension]

- Members:
    - name:      instance name
    - xform:      transform
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
    std::vector<mtl_material> materials;
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


### Struct obj_exception

~~~ .cpp
struct obj_exception : std::exception {
    obj_exception(const std::string& errmsg); 
    virtual const char* what() const throw(); 
}
~~~

IO Exception.

- Members:
    - obj_exception():      constructor with error message
    - what():      retieval of error message


### Function load_obj()

~~~ .cpp
obj* load_obj(const std::string& filename, bool flip_texcoord = true);
~~~

Load OBJ

- Parameters:
    - filename: filename
    - flip_texcoord: whether to flip the v coordinate
- Return:
    - obj

### Function load_mtl()

~~~ .cpp
std::vector<mtl_material> load_mtl(const std::string& filename);
~~~

Load MTL

- Parameters:
    - filename: filename
- Return:
    - loaded materials

### Function save_obj()

~~~ .cpp
void save_obj(
    const std::string& filename, const obj* asset, bool flip_texcoord = true);
~~~

Save OBJ

- Parameters:
    - filename: filename
    - asset: obj data to save
    - flip_texcoord: whether to flip the v coordinate

### Function save_mtl()

~~~ .cpp
void save_mtl(
    const std::string& filename, const std::vector<mtl_material>& materials);
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

