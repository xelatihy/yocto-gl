# Yocto/Shape: Shape utilities

Yocto/Shape is a collection of utilities for manipulating shapes in 3D
graphics, with a focus on triangle and quad meshes.
Yocto/Shape is implemented in `yocto_shape.h` and `yocto_shape.cpp`.

## Shape representation

Yocto/Shape supports shapes defined as collection of either points, lines,
triangles and quads. All functions have overrides for all element types
if appropriate.

Shapes are represented as indexed meshes, with arbitrary properties for
each vertex. Each vertex property is stored as a separate array,
and shape elements are stored as arrays of indices to faces.
For element parametrization, we follow [Yocto/Geometry](yocto_geometry.md).
Most functions take vertex properties and elements explicitly, without relying
on a mesh data structure.

Vertex data is stored as `vector<vecXf>`, while element indices are stored
as `vector<vec3i>`, `vector<vec4i>`, `vector<vec2i>`, `vector<int>`
for triangle meshes, quad meshes, line stets and point sets respectively.

```cpp
auto triangles = vector<vec3i>{...};   // triangle indices
auto positions = vector<vec3f>{...};   // vertex positions
auto texcoords = vector<vec2f>{...};   // vertex uvs
```

Face-varying shapes are also supported by specifying separate element indices
for each vertex property, with arrays of vertex properties possibly of different
length. This makes sure that any topology can be represented. For now, only
face-varying quads are supported.

```cpp
auto quadspos = vector<vec4i>{...};        // quads indices for positions
auto positions = vector<vec3f>{...};       // vertex positions
auto quadstexcoords = vector<vec4i>{...};  // quads indices for uvs
auto texcoords = vector<vec2f>{...};       // vertex uvs
```

## Vertex properties

Yocto/Shape provides many facilities to compute vertex properties for indexed
elements. Use `compute_normals(...)` to compute vertex normals for triangle and
quad meshes, and `compute_tangents(...)` for line tangents.
Use `compute_skinning(...)` to apply linear-blend skinning.
Use `compute_tangents_spaces(...)` to compute tangents spaces for each ech
meshes.

```cpp
auto triangles = vector<vec3i>{...};   // triangle indices
auto positions = vector<vec3f>{...};   // vertex positions
auto texcoords = vector<vec2f>{...};   // vertex uvs

auto normals = compute_normals(triangles,positions);   // vertex normals
auto tangsp = compute_tangent_spaces(triangles, positions, normals, texcoords);

auto weights = vector<vec4f>{...};   // skinning weights for 4 bones per vertex
auto joints  = vector<vec4i>{...};   // bine indices for 4 bones per vertex
auto frames  = vector<frame3f>{...}; // bone frames
auto [skinned_pos, skinned_norm] = compute_skinning(positions, normals,
   weights, joints, frames);       // skinned positions ans normals
```

## Flipping and aligning

Yocto/Shape provides functions to correct shapes that have inconsistent
orientations or normals. Use `flip_normals(normals)` to flip all mesh normals.
Use `flip_triangles(triangles)` and `flip_quads(quads)` to change face
orientations. Use `align_vertices(positions,alignment)` to align vertex
positions to the main axes.

```cpp
auto triangles = vector<vec3i>{...};   // triangle indices
auto positions = vector<vec3f>{...};   // vertex positions
auto positions = vector<vec3f>{...};   // vertex normals

triangles = flip_triangles(triangles); // flip faces
normals = flip_normals(normals);       // flip normals

// align positions to the origin along the y axis only
positions = align_vertices(positions, {0,1,0});
```

## Edges and adjacencies

Use `get_edges(triangles)` amd `get_edges(quads)` to get a list of unique edges
for a triangle or quads mesh.

```cpp
auto triangles = vector<vec3i>{...};  // triangle indices
auto edges = get_edges(triangles);    // edge indices
```

Internally, these functions use an `edge_map`, that is a dictionary that has
pairs of vertex ids as keys and an edge index as value.
Two opposing half-edges have the same representation in an `edge_map`,
making it useful in tesselation algorithms to avoid cracks.
In Yocto/Shape, edge maps also stores the number of incident faces per edge,
so that we can determine which edges belong to the boundary.

```cpp
auto triangles = vector<vec3i>{...};  // triangle indices
auto emap = make_edge_map(triangles); // edge map
auto edges = get_edges(emap);         // edge indices
for(auto& edge : edges)               // iterate over edges
   print(edge_index(emap, edge));     // get edge indices
auto boundary = get_boundary(emap);   // get unsorted boundary edges
```

## Ray-intersection and point-overlap

Yocto/Shape provides ray-scene intersection for points, lines, triangles and
quads accelerated by a BVH data structure. Our BVH is written for
minimal code and not maximum speed, but still gives fast-enough results.
See [Yocto/Geometry](yocto_geometry.md) for intersection parametrization.

The BVH tree is stored in a `bvh_tree` struct. The tree stored an array of
nodes and an array of element indices. Each node in the tree has references
to either other nodes or elements. References are represented as indices
in the nodes or elements arrays. Nodes indices refer to the nodes array,
for internal nodes, or the element arrays, for leaf nodes.
The BVH does not store shape data, which is instead passed explicitly
to all calls.

BVH nodes contain their bounds, indices to the BVH arrays of either
primitives or internal nodes, node element type,
and the split axis. Leaf and internal nodes are identical, except that
indices refer to primitives for leaf nodes or other nodes for internal nodes.

The BVH is initialized with `make_triangles_bvh(bvh,triangles,positions)` for
triangles, `make_quads_bvh(bvh,quads,positions)` for quads,
`make_lines_bvh(bvh,lines,positions,radius)` for lines, and
`make_points_bvh(bvh,points,positions,radius)` for points.

```cpp
auto triangles = vector<vec3i>{};              // mesh data
auto positions = vector<vec3f>{};
auto bvh = bvh_tree{};
make_triangles_bvh(bvh, triangles, positions;  // BVH construction
```

Intersect and overlap functions return a `bvh_intersection` that bundles
the intersection distance, the intersected element index and uvs,
and a `hit` flag that signals whether an element was hit.

`intersect_XXX_bvh(...)` computed intersections between rays and shapes.
Use `intersect_triangles_bvh(bvh,triangles,positions, ray)` for
triangles, `intersect_quads_bvh(bvh,quads,positions)` for quads,
`intersect_lines_bvh(bvh,lines,positions,radius,ray)` for lines, and
`intersect_points_bvh(bvh,points,positions,radius,ray)` for points.

```cpp
auto ray = ray3f{...};
auto isec = intersect_triangles(bvh, triangles, positions, ray);// intersection
if(isec.hit) print_info(isec.element, isec.uv, isec.distance);
else print_info("no hit");
```

`overlap_XXX_bvh(...)` checks whether a shape overlaps a point within a
given maximum distance and returns the distance, element and uv of the
closest element.
Use `overlap_triangles_bvh(bvh,triangles,positions, ray)` for
triangles, `overlap_quads_bvh(bvh,quads,positions)` for quads,
`overlap_lines_bvh(bvh,lines,positions,radius,ray)` for lines, and
`overlap_points_bvh(bvh,points,positions,radius,ray)` for points.

```cpp
// point overlaps
auto pt = vec3f{...}; auto max_dist = float{...};
auto ovr = overlap_triangles(bvh, triangles, positions, pt, mat_dist);// overlap
if(ovr.hit) print_info(ovrl.element, ovrl.uv, ovrl.distance);
else print_info("no overlap");
```

If vertices have moved little, BVHs can be updated instead of fully rebuild.
Use `update_triangles_bvh(bvh,triangles,positions)` for
triangles, `update_quads_bvh(bvh,quads,positions)` for quads,
`update_lines_bvh(bvh,lines,positions,radius)` for lines, and
`update_points_bvh(bvh,points,positions,radius)` for points.

```cpp
positions[...] = {...};                           // update positions
update_triangles_bvh(bvh, triangles, positions);  // update BVH
```

## Nearest neighbors using a hash grid

Nearest neighbors queries are computed by building a sparse hash grid
defined as `hash_grid`. The grid is created by specifying a cell size for the
underlying volumetric grid. Each cell stores the list of point indices that
are present in that cell. To save memory, the grid is represented sparsely,
using a dictionary, so that only cells with at least one vertex are defined.

Initialize a hash grid with `make_hash_grid(positions,size)`.
Use `find_neighbors(grid, neighbors, position, max_radius)` to find
nearest neighbors.

```cpp
auto positions = vector<vec3f>{...};               // point positions
auto grid = make_hash_grid(positions,cell_size);   // create hash grid

auto pt = vec3f{...}; auto max_dist = float{...};  // query point and dist
auto neighbors = vector<int>{};                    // neighbor buffer
find_neighbors(grid, neighbors, pt, max_dist);     // find neighbors by pos

find_neighbors(grid, neighbors, id, max_dist);     // find neighbors by id
```

## Element conversions and grouping

Yocto/Shape support conversion between shape elements.
Use `quads_to_triangles(quads)` to convert quads to triangles and
`triangles_to_quads(triangles)` to convert triangles to degenerate quads.
Use `bezier_to_lines(lines)` to convert Bézier segments to lines using three
lines for each Bézier segment.

```cpp
auto quads = vector<vec4i>{...};
auto triangles = quads_to_triangles(quads);  // convert quads to triangles
```

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm, texcoord and colors.
std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
split_facevarying(const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
const vector<vec3f>& normals, const vector<vec2f>& texcoords);

// Split primitives per id
vector<vector<vec2i>> ungroup_lines(
const vector<vec2i>& lines, const vector<int>& ids);
vector<vector<vec3i>> ungroup_triangles(
const vector<vec3i>& triangles, const vector<int>& ids);
vector<vector<vec4i>> ungroup_quads(
const vector<vec4i>& quads, const vector<int>& ids);

// Weld vertices within a threshold.
std::pair<vector<vec3f>, vector<int>> weld_vertices(
const vector<vec3f>& positions, float threshold);
std::pair<vector<vec3i>, vector<vec3f>> weld_triangles(
const vector<vec3i>& triangles, const vector<vec3f>& positions,
float threshold);
std::pair<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
const vector<vec3f>& positions, float threshold);

// Merge shape elements
void merge_lines(
vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts);
void merge_triangles(vector<vec3i>& triangles,
const vector<vec2i>& merge_triangles, int num_verts);
void merge_quads(
vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts);
void merge_lines(vector<vec2i>& lines, vector<vec3f>& positions,
vector<vec3f>& tangents, vector<vec2f>& texcoords, vector<float>& radius,
const vector<vec2i>& merge_lines, const vector<vec3f>& merge_positions,
const vector<vec3f>& merge_tangents,
const vector<vec2f>& merge_texturecoords,
const vector<float>& merge_radius);
void merge_triangles(vector<vec3i>& triangles, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vector<vec2i>& merge_triangles, const vector<vec3f>& merge_positions,
const vector<vec3f>& merge_normals,
const vector<vec2f>& merge_texturecoords);
void merge_quads(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vector<vec4i>& merge_quads, const vector<vec3f>& merge_positions,
const vector<vec3f>& merge_normals,
const vector<vec2f>& merge_texturecoords);

// Merge quads and triangles
void merge_triangles_and_quads(
vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles);

## Shape subdivision

// Subdivide lines by splitting each line in half.
std::pair<vector<vec2i>, vector<float>> subdivide_lines(
const vector<vec2i>& lines, const vector<float>& vert, int level);
std::pair<vector<vec2i>, vector<vec2f>> subdivide_lines(
const vector<vec2i>& lines, const vector<vec2f>& vert, int level);
std::pair<vector<vec2i>, vector<vec3f>> subdivide_lines(
const vector<vec2i>& lines, const vector<vec3f>& vert, int level);
std::pair<vector<vec2i>, vector<vec4f>> subdivide_lines(
const vector<vec2i>& lines, const vector<vec4f>& vert, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
std::pair<vector<vec3i>, vector<float>> subdivide_triangles(
const vector<vec3i>& triangles, const vector<float>& vert, int level);
std::pair<vector<vec3i>, vector<vec2f>> subdivide_triangles(
const vector<vec3i>& triangles, const vector<vec2f>& vert, int level);
std::pair<vector<vec3i>, vector<vec3f>> subdivide_triangles(
const vector<vec3i>& triangles, const vector<vec3f>& vert, int level);
std::pair<vector<vec3i>, vector<vec4f>> subdivide_triangles(
const vector<vec3i>& triangles, const vector<vec4f>& vert, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
std::pair<vector<vec4i>, vector<float>> subdivide_quads(
const vector<vec4i>& quads, const vector<float>& vert, int level);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_quads(
const vector<vec4i>& quads, const vector<vec2f>& vert, int level);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_quads(
const vector<vec4i>& quads, const vector<vec3f>& vert, int level);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_quads(
const vector<vec4i>& quads, const vector<vec4f>& vert, int level);
// Subdivide beziers by splitting each segment in two.
std::pair<vector<vec4i>, vector<float>> subdivide_beziers(
const vector<vec4i>& beziers, const vector<float>& vert, int level);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_beziers(
const vector<vec4i>& beziers, const vector<vec2f>& vert, int level);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_beziers(
const vector<vec4i>& beziers, const vector<vec3f>& vert, int level);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_beziers(
const vector<vec4i>& beziers, const vector<vec4f>& vert, int level);
// Subdivide quads using Carmull-Clark subdivision rules.
std::pair<vector<vec4i>, vector<float>> subdivide_catmullclark(
const vector<vec4i>& quads, const vector<float>& vert, int level,
bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
const vector<vec4i>& quads, const vector<vec2f>& vert, int level,
bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
const vector<vec4i>& quads, const vector<vec3f>& vert, int level,
bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
const vector<vec4i>& quads, const vector<vec4f>& vert, int level,
bool lock_boundary = false);

## Shape sampling

// Pick a point in a point set uniformly.
int sample_points(int npoints, float re);
int sample_points(const vector<float>& cdf, float re);
vector<float> sample_points_cdf(int npoints);

// Pick a point on lines uniformly.
std::pair<int, float> sample_lines(
const vector<float>& cdf, float re, float ru);
vector<float> sample_lines_cdf(
const vector<vec2i>& lines, const vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
std::pair<int, vec2f> sample_triangles(
const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_triangles_cdf(
const vector<vec3i>& triangles, const vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
std::pair<int, vec2f> sample_quads(
const vector<float>& cdf, float re, const vec2f& ruv);
std::pair<int, vec2f> sample_quads(const vector<vec4i>& quads,
const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_quads_cdf(
const vector<vec4i>& quads, const vector<vec3f>& positions);

// Samples a set of points over a triangle/quad mesh uniformly. Returns pos,
// norm and texcoord of the sampled points.
void sample_triangles(vector<vec3f>& sampled_positions,
vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
const vector<vec3i>& triangles, const vector<vec3f>& positions,
const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
int seed = 7);
void sample_quads(vector<vec3f>& sampled_positions,
vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
const vector<vec4i>& quads, const vector<vec3f>& positions,
const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
int seed = 7);

## Shape IO

// Load/save a shape as indexed meshes
[[nodiscard]] bool load_shape(const string& filename, vector<int>& points,
vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
vector<vec3f>& colors, vector<float>& radius, string& error,
bool flip_texcoords = true);
[[nodiscard]] bool save_shape(const string& filename, const vector<int>& points,
const vector<vec2i>& lines, const vector<vec3i>& triangles,
const vector<vec4i>& quads, const vector<vec3f>& positions,
const vector<vec3f>& normals, const vector<vec2f>& texcoords,
const vector<vec3f>& colors, const vector<float>& radius, string& error,
bool ascii = false, bool flip_texcoords = true);

// Load/save a facevarying shape
[[nodiscard]] bool load_fvshape(const string& filename, vector<vec4i>& quadspos,
vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_fvshape(const string& filename,
const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
const vector<vec3f>& normals, const vector<vec2f>& texcoords, string& error,
bool ascii = false, bool flip_texcoords = true);

## Shape stats

// Get mesh statistics for printing
vector<string> shape_stats(const vector<int>& points,
const vector<vec2i>& lines, const vector<vec3i>& triangles,
const vector<vec4i>& quads, const vector<vec4i>& quadspos,
const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
const vector<vec3f>& positions, const vector<vec3f>& normals,
const vector<vec2f>& texcoords, const vector<vec3f>& colors,
const vector<float>& radius, bool verbose = false);

## Procedural shapes

// Make a plane.
void make_rect(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1});
void make_bulged_rect(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1}, float radius = 0.3);
void make_recty(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1});
// Make a box.
void make_box(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
const vec3f& uvscale = {1, 1, 1});
void make_rounded_box(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
const vec3f& uvscale = {1, 1, 1}, float radius = 0.3);
// Make a quad stack
void make_rect_stack(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
const vec3f& uvscale = {1, 1, 1});
// Make a floor.
void make_floor(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
const vec2f& uvscale = {10, 10});
void make_bent_floor(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
const vec2f& uvscale = {10, 10}, float bent = 0.5);
// Make a sphere.
void make_sphere(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
float scale = 1, float uvscale = 1);
// Make a sphere.
void make_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {32, 32}, float scale = 1,
const vec2f& uvscale = {1, 1});
// Make a sphere with slipped caps.
void make_capped_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {32, 32}, float scale = 1,
const vec2f& uvscale = {1, 1}, float height = 0.3);
// Make a disk
void make_disk(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
float scale = 1, float uvscale = 1);
// Make a bulged disk
void make_bulged_disk(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
float scale = 1, float uvscale = 1, float height = 0.3);
// Make a uv disk
void make_uvdisk(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {32, 32}, float scale = 1,
const vec2f& uvscale = {1, 1});
// Make a uv cylinder
void make_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1},
const vec3f& uvscale = {1, 1, 1});
// Make a rounded uv cylinder
void make_rounded_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1},
const vec3f& uvscale = {1, 1, 1}, float radius = 0.3);
// Make a plane in the xz plane.
void make_yrect(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1});
void make_bulged_yrect(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1}, float radius = 0.3);

// Make a facevarying rect
void make_fvrect(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
const vec2f& uvscale = {1, 1});
// Make a facevarying box
void make_fvbox(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords,
const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
const vec3f& uvscale = {1, 1, 1});
void make_fvsphere(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
float scale = 1, float uvscale = 1);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
void make_lines(vector<vec2i>& lines, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
int num = 65536, const vec2i& steps = {4, 65536},
const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1},
const vec2f& rad = {0.001, 0.001});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
void make_point(vector<int>& points, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
float point_radius);
void make_points(vector<int>& points, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
int num = 65536, float uvscale = 1, float point_radius = 0.001);
void make_random_points(vector<int>& points, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
int num = 65536, const vec3f& size = {1, 1, 1}, float uvscale = 1,
float point_radius = 0.001, uint64_t seed = 17);

// Predefined meshes
void make_monkey(
vector<vec4i>& quads, vector<vec3f>& positions, float scale = 1);
void make_quad(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_quady(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_cube(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_fvcube(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_geosphere(
vector<vec3i>& triangles, vector<vec3f>& positions, float scale = 1);

// Make a hair ball around a shape. Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (strength/number)
// rotation: rotation added to hair (angle/strength)
void make_hair(vector<vec2i>& lines, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
const vector<vec3i>& striangles, const vector<vec4i>& squads,
const vector<vec3f>& spos, const vector<vec3f>& snorm,
const vector<vec2f>& stexcoord, const vec2i& steps = {8, 65536},
const vec2f& length = {0.1, 0.1}, const vec2f& rad = {0.001, 0.001},
const vec2f& noise = {0, 10}, const vec2f& clump = {0, 128},
const vec2f& rotation = {0, 0}, int seed = 7);

// Thickens a shape by copying the shape content, rescaling it and flipping its
// normals. Note that this is very much not robust and only useful for trivial
// cases.
void make_shell(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, float thickness);

// Make a heightfield mesh.
void make_heightfield(vector<vec4i>& quads, vector<vec3f>& positions,
vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& size,
const vector<float>& height);
