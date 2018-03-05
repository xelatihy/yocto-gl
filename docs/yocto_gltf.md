# Yocto/glTF: High-level glTF support for Yocto/GL

Extension of Yocto/GL that provides a complete high-level interface for
Khronos glTF format. Uses the low-level parser from Yocto/GL.
The high-level interface that is more useful if an aplication needs direct
access to shapes, textures and animations and does not want to deal with
all the intricacies of the format.

Known limitations are: (1) skinning matrices are always in world space
(waiting for the spec to be updated); (2) spline based animation is not
fully implemented yet (waiting for official demo models).

This library depends on `yocto_gl.h` and its dependencies.


## Usage

1. load a group of scens with `load_scenes()`
2. look at the `scene` data structures for access to individual elements
3. to support animation, use `update_animated_transforms()`
4. to support skinning, use `get_skin_transforms()`
5. for morphing, use `compute_morphing_deformation()`
6. can also manipulate the scene by adding missing data with `add_XXX()`
   functions
7. for rendering scenes, use `get_scene_cameras()` and
   `get_scene_instances()` that avoid the need for explicirtly walking
   the glTF node hierarchy
8. use `save_scenes()` to write the data to disk
9. use `convert_to_specgloss()` to convert materials to spec-gloss


## History

- v 0.1.0: initial release after refactoring

## API Documentation

#### Struct gltf_camera

~~~ .cpp
struct gltf_camera {
    std::string name = "";
    bool ortho = false;
    float aspect = 1;
    float yfov = 2 * std::atan(0.5f);
    float near = 0;
    float far = 0;
    float focus = 1;
    float aperture = 0;
}
~~~

Camera

- Members:
    - name:      name
    - ortho:      orthographic
    - aspect:      aspect ratio
    - yfov:      vertical fov (perspective) or size (orthographic)
    - near:      near plane (0 for default)
    - far:      far plane (0 for default)
    - focus:      focus distance (extension not implemented yet)
    - aperture:      lens aperture (extension not implemented yet)


#### Struct gltf_texture

~~~ .cpp
struct gltf_texture {
    std::string name = "";
    std::string path = "";
    image4b ldr;
    image4f hdr;
    int width() const; 
    int height() const; 
}
~~~

Texture

- Members:
    - name:      name
    - path:      path
    - ldr:      8-bit data
    - hdr:      float data
    - width():      get texture width
    - height():      get texture height


#### Enum gltf_texture_wrap

~~~ .cpp
enum struct gltf_texture_wrap {
    repeat = 1,
    clamp = 2,
    mirror = 3,
}
~~~

Texture wrap mode

- Values:
    - repeat:      repeat
    - clamp:      clamp
    - mirror:      mirror


#### Enum gltf_texture_filter

~~~ .cpp
enum struct gltf_texture_filter {
    linear = 1,
    nearest = 2,
    linear_mipmap_linear = 3,
    nearest_mipmap_nearest = 4,
    linear_mipmap_nearest = 5,
    nearest_mipmap_linear = 6,
}
~~~

Texture filter mode

- Values:
    - linear:      linear
    - nearest:      nearest
    - linear_mipmap_linear:      linear mipmap linear
    - nearest_mipmap_nearest:      nearest mipmap nearest
    - linear_mipmap_nearest:      linear mipmap nearest
    - nearest_mipmap_linear:      nearest mipmap linear


#### Struct gltf_texture_info

~~~ .cpp
struct gltf_texture_info {
    gltf_texture_wrap wrap_s = gltf_texture_wrap::repeat;
    gltf_texture_wrap wrap_t = gltf_texture_wrap::repeat;
    gltf_texture_filter filter_mag = gltf_texture_filter::linear;
    gltf_texture_filter filter_min = gltf_texture_filter::linear_mipmap_linear;
    float scale = 1;
    bool is_default() const; 
}
~~~

Texture information

- Members:
    - wrap_s:      wrap mode for s coordinate
    - wrap_t:      wrap mdoe for t coordinate
    - filter_mag:      filter mode
    - filter_min:      filter mode
    - scale:      texture strength (occlusion and normal)
    - is_default():      check if it is default


#### Struct gltf_material_metallic_roughness

~~~ .cpp
struct gltf_material_metallic_roughness {
    vec3f base = {0, 0, 0};
    float opacity = 1;
    float metallic = 0;
    float roughness = 0;
    gltf_texture* base_txt = nullptr;
    gltf_texture* metallic_txt = nullptr;
    gltf_texture_info* base_txt_info = nullptr;
    gltf_texture_info* metallic_txt_info = nullptr;
}
~~~

Material PBR metallic roughness

- Members:
    - base:      base color
    - opacity:      opacity
    - metallic:      metallic factor
    - roughness:      metallic roughness
    - base_txt:      base texture (kb.x, kb.y, kb.z, op)
    - metallic_txt:      metallic-roughness texture (n/a, roughness, metallic, n/a)
    - base_txt_info:      texture information for base_txt
    - metallic_txt_info:      texture information for metallic_txt


#### Struct gltf_material_specular_glossiness

~~~ .cpp
struct gltf_material_specular_glossiness {
    vec3f diffuse = {0, 0, 0};
    float opacity = 1;
    vec3f specular = {0, 0, 0};
    float glossiness = 1;
    gltf_texture* diffuse_txt = nullptr;
    gltf_texture* specular_txt = nullptr;
    gltf_texture_info* diffuse_txt_info = nullptr;
    gltf_texture_info* specular_txt_info = nullptr;
}
~~~

Material PBR specular glossiness

- Members:
    - diffuse:      diffuse color and opacity
    - opacity:      opacity
    - specular:      specular color (spec.x, spec.y, spec.z, opacity)
    - glossiness:      specular glossiness
    - diffuse_txt:      diffuse texture (diff.x, diff.y, diff.z, opacity)
    - specular_txt:      specular-glossiness texture (spec.x, spec.y, spec.z, gloss)
    - diffuse_txt_info:      texture information for base_txt
    - specular_txt_info:      texture information for metallic_txt


#### Struct gltf_material

~~~ .cpp
struct gltf_material {
    std::string name = "";
    vec3f emission = {0, 0, 0};
    gltf_texture* emission_txt = nullptr;
    gltf_texture_info* emission_txt_info = nullptr;
    gltf_material_metallic_roughness* metallic_roughness = nullptr;
    gltf_material_specular_glossiness* specular_glossiness = nullptr;
    gltf_texture* occlusion_txt = nullptr;
    gltf_texture* normal_txt = nullptr;
    gltf_texture_info* occlusion_txt_info = nullptr;
    gltf_texture_info* normal_txt_info = nullptr;
    bool double_sided = true;
}
~~~

Material

glTF 2.0 has two physically-based aterial models: pbrMetallicRoughness
and pbrSpecularGlossiness, the latter as an extension. Here we support both
by including which one is defined. While it would have been more appropriate
to convert them, this requires a rewrite of texture data which w prefer to
avoid.

- Members:
    - name:      name
    - emission:      emission color
    - emission_txt:      emissive texture reference
    - emission_txt_info:      texture information for normal_txt
    - metallic_roughness:      metallic roughnesss
    - specular_glossiness:      specular glossiness
    - occlusion_txt:      occlusion texture
    - normal_txt:      normal texture
    - occlusion_txt_info:      texture information for collusion_txt
    - normal_txt_info:      texture information for normal_txt
    - double_sided:      double sided


#### Struct gltf_shape_morph

~~~ .cpp
struct gltf_shape_morph {
    std::vector<vec3f> pos;
    std::vector<vec3f> norm;
    std::vector<vec3f> tangsp;
    float weight = 0;
}
~~~

Morph information for shapes

- Members:
    - pos:      morph position
    - norm:      morph normal
    - tangsp:      morph tangent
    - weight:      default weight (the same for each shape in a mesh)


#### Struct gltf_shape

~~~ .cpp

struct gltf_shape {
    std::string name = "";
    gltf_material* mat = nullptr;
    std::vector<vec3f> pos;
    std::vector<vec3f> norm;
    std::vector<vec2f> texcoord;
    std::vector<vec2f> texcoord1;
    std::vector<vec4f> color;
    std::vector<float> radius;
    std::vector<vec4f> tangsp;
    std::vector<vec4f> skin_weights;
    std::vector<vec4i> skin_joints;
    std::vector<int> points;
    std::vector<vec2i> lines;
    std::vector<vec3i> triangles;
    std::vector<gltf_shape_morph*> morph_targets;
}
~~~

Primitives

- Members:
    - name:      name of the mesh that enclosed it
    - mat:      material reference
    - pos:      vertex position
    - norm:      vertex normal
    - texcoord:      vertex texcoord
    - texcoord1:      vertex additional texcoord
    - color:      vertex color
    - radius:      vertex radius
    - tangsp:      vertex tangent space
    - skin_weights:      vertex skinning weights
    - skin_joints:      vertex skinning joint indices
    - points:      point elements
    - lines:      line elements
    - triangles:      triangle elements
    - morph_targets:      morph targets


#### Struct gltf_mesh

~~~ .cpp
struct gltf_mesh {
    std::string name = "";
    std::string path = "";
    std::vector<gltf_shape*> shapes;
}
~~~

Gltf mesh.

- Members:
    - name:      name
    - path:      path (only used when writing files on disk with glTF)
    - shapes:      primitives references


#### Struct gltf_node

~~~ .cpp
struct gltf_node {
    std::string name = "";
    gltf_camera* cam = nullptr;
    gltf_mesh* msh = nullptr;
    gltf_skin* skn = nullptr;
    std::vector<gltf_node*> children;
    mat4f matrix = identity_mat4f;
    quat4f rotation = {0, 0, 0, 1};
    vec3f scale = {1, 1, 1};
    vec3f translation = {0, 0, 0};
    std::vector<float> morph_weights;
    gltf_node* parent = nullptr;
    mat4f xform() const; 
    mat4f local_xform() const; 
    mat4f skin_xform() const; 
    mat4f _xform = identity_mat4f;
    mat4f _local_xform = identity_mat4f;
    mat4f _skin_xform = identity_mat4f;
}
~~~

Node in the hierarchy.

- Members:
    - name:      name
    - cam:      camera reference
    - msh:      mesh reference
    - skn:      mesh reference
    - children:      children
    - matrix:      A floating-point 4x4 transformation matrix stored in column-major order.
    - rotation:      The node's unit quaternion rotation in the order (x, y, z, w), where w
     is the scalar.
    - scale:      The node's non-uniform scale.
    - translation:      The node's translation.
    - morph_weights:      morph target weights
    - parent:      parent node (computed during update_node_hierarchy())
    - xform():      transform (computed during update_transforms())
    - local_xform():      local transform (computed during update_transforms())
    - skin_xform():      skin transform (computed during update_transforms())
    - _xform:      transform (computed during update_transforms())
    - _local_xform:      local transform (computed during update_transforms())
    - _skin_xform:      skin transform (computed during update_transforms())


#### Enum gltf_animation_interpolation

~~~ .cpp
enum struct gltf_animation_interpolation {
    linear = 0,
    step = 1,
    catmull_rom = 2,
    cubic = 3,
}
~~~

Animation Interpolation

- Values:
    - linear:      linear
    - step:      step function
    - catmull_rom:      catmull-rom spline
    - cubic:      cubic bezier spline


#### Struct gltf_animation

~~~ .cpp
struct gltf_animation {
    gltf_animation_interpolation interp = gltf_animation_interpolation::step;
    std::vector<gltf_node*> nodes;
    std::vector<float> time;
    std::vector<vec3f> translation;
    std::vector<quat4f> rotation;
    std::vector<vec3f> scale;
    std::vector<std::vector<float>> morph_weights;
}
~~~

Keyframe data.

- Members:
    - interp:      Interpolation
    - nodes:      Target nodes
    - time:      Times
    - translation:      Translation
    - rotation:      Rotation
    - scale:      Scale
    - morph_weights:      Weights for morphing


#### Struct gltf_animation_group

~~~ .cpp
struct gltf_animation_group {
    std::string name;
    std::string path = "";
    std::vector<gltf_animation*> animations;
}
~~~

Animation

- Members:
    - name:      Name
    - path:      path (only used when writing files on disk with glTF)
    - animations:      Times


#### Struct gltf_skin

~~~ .cpp
struct gltf_skin {
    std::string name = "";
    std::string path = "";
    std::vector<mat4f> pose_matrices;
    std::vector<gltf_node*> joints;
    gltf_node* root = nullptr;
}
~~~

Skin

- Members:
    - name:      name
    - path:      path (only used when writing files on disk with glTF)
    - pose_matrices:      inverse bind matrix
    - joints:      joints
    - root:      skeleton root node


#### Struct gltf_scene

~~~ .cpp
struct gltf_scene {
    std::string name = "";
    std::vector<gltf_node*> nodes;
}
~~~

Gltf scene

- Members:
    - name:      name
    - nodes:      instances


#### Struct gltf_scene_group

~~~ .cpp
struct gltf_scene_group {
    gltf_scene* default_scene = nullptr;
    std::vector<gltf_camera*> cameras;
    std::vector<gltf_material*> materials;
    std::vector<gltf_texture*> textures;
    std::vector<gltf_mesh*> meshes;
    std::vector<gltf_scene*> scenes;
    std::vector<gltf_node*> nodes;
    std::vector<gltf_animation_group*> animations;
    std::vector<gltf_skin*> skins;
}
~~~

Gltf model. Objects are shared between scenes.
Scenes and nodes are missing for mesh-only assets.

- Members:
    - default_scene:      default scene (null if not present)
    - cameras:      cameras
    - materials:      materials
    - textures:      textures
    - meshes:      meshes
    - scenes:      scenes
    - nodes:      nodes
    - animations:      nodes
    - skins:      skins


#### Function load_scenes()

~~~ .cpp
gltf_scene_group* load_scenes(
    const std::string& filename, bool load_textures, bool skip_missing = true);
~~~

Load scene

- Parameters:
    - filename: filename
    - load_textures: whether to load textures (default to false)
    - skip_missing: whether to skip missing buffers and textures
- Returns:
    - scene (nullptr on error)

#### Function save_scenes()

~~~ .cpp
void save_scenes(const std::string& filename, const std::string& buffer_uri,
    const gltf_scene_group* scn, bool save_textures,
    bool separate_buffers = false);
~~~

Save scene

- Parameters:
    - filename: filename
    - buffer_uri: name of the main buffer
    - scn: scene data to save
    - save_textures: whether to save textures (default to false)
    - separate_buffers: save separate buffers for each mesh

#### Function update_node_hierarchy()

~~~ .cpp
void update_node_hierarchy(gltf_scene_group* scn);
~~~

Update node hierarchy

#### Function update_transforms()

~~~ .cpp
void update_transforms(gltf_scene_group* scn);
~~~

Update node trasforms

#### Function update_animated_transforms()

~~~ .cpp
void update_animated_transforms(gltf_scene_group* scns, float time);
~~~

Update animated node

#### Function get_mesh_nodes()

~~~ .cpp
std::vector<gltf_node*> get_mesh_nodes(const gltf_scene* scn);
~~~

Get a list of nodes with meshes

#### Function get_camera_nodes()

~~~ .cpp
std::vector<gltf_node*> get_camera_nodes(const gltf_scene* scn);
~~~

Get a list of nodes with cameras

#### Function get_animation_bounds()

~~~ .cpp
vec2f get_animation_bounds(const gltf_scene_group* scn);
~~~

Animation times

#### Function get_skin_transforms()

~~~ .cpp
std::vector<mat4f> get_skin_transforms(const gltf_skin* sk, const mat4f& xform);
~~~

Skin transforms (local-to-object) from the node transform that instances the
skin

#### Function compute_morphing_deformation()

~~~ .cpp
void compute_morphing_deformation(const gltf_shape* shp,
    const std::vector<float>& weights, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec4f>& tangsp);
~~~

Compute shape morphing

#### Function compute_scene_bounds()

~~~ .cpp
bbox3f compute_scene_bounds(const gltf_scene_group* scn);
~~~

Computes a scene bounding box

#### Function add_normals()

~~~ .cpp
void add_normals(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_radius()

~~~ .cpp
void add_radius(gltf_scene_group* scn, float radius);
~~~

Add missing data to the scene.

#### Function add_tangent_space()

~~~ .cpp
void add_tangent_space(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_nodes()

~~~ .cpp
void add_nodes(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_scene()

~~~ .cpp
void add_scene(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_texture_data()

~~~ .cpp
void add_texture_data(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_names()

~~~ .cpp
void add_names(gltf_scene_group* scn);
~~~

Add missing data to the scene.

#### Function add_default_cameras()

~~~ .cpp
void add_default_cameras(gltf_scene_group* scn);
~~~

Add a default camera that views the entire scene.

#### Function add_unique_path_names()

~~~ .cpp
void add_unique_path_names(
    gltf_scene_group* scns, const std::string& buffer_uri);
~~~

Set unique path names for outputting separate buffers

#### Function add_spec_gloss()

~~~ .cpp
void add_spec_gloss(gltf_scene_group* scns);
~~~

Convert materials to spec gloss

#### Function gltf_to_scenes()

~~~ .cpp
gltf_scene_group* gltf_to_scenes(const glTF* gltf, int scene_idx = -1);
~~~

Convert a gltf asset to flattened group of scene.

#### Function scenes_to_gltf()

~~~ .cpp
glTF* scenes_to_gltf(const gltf_scene_group* fl_gltf,
    const std::string& buffer_uri, bool separate_buffers = false);
~~~

Convert a flattened group of scene into a gltf. If separate_buffers,
creates a separate buffer for each each and animation and
prepend buffer_uri to its name.

#### Function validate_gltf()

~~~ .cpp
std::vector<std::pair<std::string, std::string>> validate_gltf(glTF* gltf);
~~~

Validate a gltf. Missing many validation as of this version.

