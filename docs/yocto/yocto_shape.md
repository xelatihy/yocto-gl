# Yocto/Shape: Shape utilities

Yocto/Shape is a collection of utilities for manipulating shapes in 3D
graphics, with a focus on triangle and quad meshes.
Yocto/Shape is implemented in `yocto_shape.h` and `yocto_shape.cpp`.

## Shape representation

Yocto/Shape supports shapes defined as collection of either points, lines,
triangles and quads. Most functions have overrides for all element types
when appropriate.

Shapes are represented as indexed meshes, with arbitrary properties for
each vertex. Each vertex property is stored as a separate array,
and shape elements are stored as arrays of indices to faces.
For element parametrization, we follow [Yocto/Geometry](yocto_geometry.md).

Vertex data is stored as `vector<vecXf>`, while element indices are stored
as `vector<vec3i>`, `vector<vec4i>`, `vector<vec2i>`, `vector<int>`
for triangle meshes, quad meshes, line sets and point sets respectively.

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

Throughout the library, functions may either take index and vertex arrays
directly as input and output, or may pack these array in structs if deemed
appropriate. This design tries to balance readability and generality, without
forcing a single convention that would not be appropriate everywhere.

If a higher level design is needed, [Yocto/Scene](yocto_scene.md) contains
the standalone types `scene_shape` and `scene_fvshape` to store indexed
and face-varying shapes respectively, a a collection of methods to work
on these type, that are essentially wrappers to the functionality in this
library.

Shape loading and saving is defined in [Yocto/SceneIO](yocto_sceneio.md).

## Vertex properties

Yocto/Shape provides many facilities to compute vertex properties for indexed
elements. Use `triangles_normals(...)` and `quads_normals(...)` to compute
vertex normals for triangle and quad meshes, and `line_tangents(...)` for
line tangents. Use `skin_vertices(...)` to apply linear-blend skinning.
Use `triangle_tangent_spaces(...)` to compute tangents spaces for each ech
meshes.

```cpp
auto triangles = vector<vec3i>{...};   // triangle indices
auto positions = vector<vec3f>{...};   // vertex positions
auto texcoords = vector<vec2f>{...};   // vertex uvs

auto normals = triangle_normals(triangles,positions);   // vertex normals
auto tangsp = triangle_tangent_spaces(triangles, positions, normals, texcoords);

auto weights = vector<vec4f>{...};   // skinning weights for 4 bones per vertex
auto joints  = vector<vec4i>{...};   // bine indices for 4 bones per vertex
auto frames  = vector<frame3f>{...}; // bone frames
auto [skinned_pos, skinned_norm] = skin_vertices(positions, normals,
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
auto normals   = vector<vec3f>{...};   // vertex normals

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
See [Yocto/Geometry](yocto_geometry.md) for intersection parametrization,
and [Yocto/Bvh](yocto_bvh.md) for a more comprehensive version.

The BVH tree is stored in a `shape_bvh` struct. The tree stored an array of
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
auto triangles = vector<vec3i>{...};                 // mesh data
auto positions = vector<vec3f>{...};
auto bvh = make_triangles_bvh(triangles, positions); // BVH construction
```

Intersect and overlap functions return a `bvh_intersection` that bundles
the intersection distance, the intersected element index and uvs,
and a `hit` flag that signals whether an element was hit.

`intersect_<element>_bvh(...)` computed intersections between rays and shapes.
Use `intersect_triangles_bvh(bvh,triangles,positions, ray)` for
triangles, `intersect_quads_bvh(bvh,quads,positions)` for quads,
`intersect_lines_bvh(bvh,lines,positions,radius,ray)` for lines, and
`intersect_points_bvh(bvh,points,positions,radius,ray)` for points.

```cpp
auto ray = ray3f{...};
// computes ray-triangles intersection
auto isec = intersect_triangles_bvh(bvh, triangles, positions, ray);
if(isec.hit) print_info(isec.element, isec.uv, isec.distance);
else print_info("no hit");
```

`overlap_<element>_bvh(...)` checks whether a shape overlaps a point within a
given maximum distance and returns the distance, element and uv of the
closest element.
Use `overlap_triangles_bvh(bvh, triangles, positions, ray)` for
triangles, `overlap_quads_bvh(bvh, quads, positions)` for quads,
`overlap_lines_bvh(bvh, lines, positions, radius, ray)` for lines, and
`overlap_points_bvh(bvh, points, positions, radius, ray)` for points.

```cpp
auto pt = vec3f{...}; auto max_dist = float{...};
// comnpute point-triangles overlap
auto ovr = overlap_triangles_bvh(bvh, triangles, positions, pt, mat_dist);
if(ovr.hit) print_info(ovrl.element, ovrl.uv, ovrl.distance);
else print_info("no overlap");
```

If vertices have moved little, BVHs can be updated instead of fully rebuild.
Use `update_triangles_bvh(bvh, triangles, positions)` for
triangles, `update_quads_bvh(bvh, quads, positions)` for quads,
`update_lines_bvh(bvh, lines, positions, radius)` for lines, and
`update_points_bvh(bvh, points, positions, radius)` for points.

```cpp
positions[...] = {...};                           // update positions
update_triangles_bvh(bvh, triangles, positions);  // update BVH
```

## Nearest neighbors

Nearest neighbors queries are computed by building a sparse hash grid
defined as `hash_grid`. The grid is created by specifying a cell size for the
underlying volumetric grid. Each cell stores the list of point indices that
are present in that cell. To save memory, the grid is represented sparsely,
using a dictionary, so that only cells with at least one vertex are defined.

Initialize a hash grid with `make_hash_grid(positions, size)`.
Use `find_neighbors(grid, neighbors, position, max_radius)` to find
nearest neighbors.

```cpp
auto positions = vector<vec3f>{...};               // point positions
auto grid = make_hash_grid(positions, cell_size);  // create hash grid
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

Face-varying meshes are stored by having different face indices for each
vertex propeerty. This way, every vertex property has its own topology.
Use `split_facevarying(...)` to convert to an indexed mesh. During conversion
vertices may be duplicated since the same topology is used for all vertex
properties.

```cpp
auto fvquadspos = vector<vec4i>{...};      // face-varying indices
auto fvquadsnorm = vector<vec4i>{...};     // arrays have some length
auto fvquadstexcoord = vector<vec4i>{...};
auto fvpositions = vector<vec3f>{...};     // face-varying vertices
auto fvnormals = vector<vec3f>{...};       // arrays may have different lengths
auto fvtexcoords = vector<vec2f>{...};
auto [quads, positions, normals, texcoords] = // convert to indexed mesh
   split_facevarying(fvquadspos, fvquadsnorm, fvquadstexcoord,
   fvpositions, fvnormals, fvtexcoords);
```

Yocto/Shape supports eliminating duplicate vertices in triangle and quad meshes.
All vertices within a threshold are merged in a greedy fashion, which works well
when duplicated vertices are near other while other vertices are further away.
Use `weld_triangles(triangles, positions, threshold)` to eliminate duplicated
triangle vertices and `weld_quads(quads, positions, threshold)` to eliminate
duplicated quad vertices. For lower-level algorithms, use
`weld_vertices(positions, threshold)` to group vertices together.

```cpp
auto triangles = vector<vec3i>{...};  // mesh data
auto positions = vector<vec3f>{...};
auto tolerance = 0.1f;
auto [mtriangles, mpositions] =       // remove duplicates
   weld_triangles(triangles, positions, tolerance);
```

Yocto/Shape supports splitting shapes that are tagged by ids. This is helpful
for example when drawing meshes that have per-face materials using renders
that do support one material per shape only. Use `ungroup_lines(lines,ids)`,
`ungroud_triangles(triangles,ids)` and `ungroup_quads(quads,ids)` for lines,
triangles and quads respectively.

```cpp
auto triangles = vector<vec3i>{...};  // tagged mesh with one id per face
auto ids       = vector<int>{...};
auto split = ungroup_triangles(triangles, ids); // returns list of meshes
```

Yocto/Shape supports merging shape elements. This is useful, for example, when
building up shapes from parts. The merged shapes are just concatenation of
the individual shape without vertex merging. Use `merge_lines(...)` for lines,
`merge_triangles(...)` for triangles and `merge_quads(...)` for quads.

```cpp
auto triangles = vector<vec3i>{...};   // initial shape
auto positions = vector<vec3f>{...};
auto normals = vector<vec3f>{...};
auto texcoords = vector<vec2f>{...};
auto mtriangles = vector<vec3i>{...};   // shape to be merged
auto mpositions = vector<vec3f>{...};
auto mnormals = vector<vec3f>{...};
auto mtexcoords = vector<vec2f>{...};
// merge mtriangles into triangles in-placee
merge_triangles(triangles, positions, normals, texcoords,
  mtriangles, mpositions, mnormals, mtexcoords);
```

You can also merge triangles and quads together in other to have one
primitive only. Use `merge_triangles_and_quads(triangles, quads, force_triangles)`
to merge elements in-place. The algorithms will output quads if present or
triangles if not unless `force_triangles` is used.

## Shape subdivision

Yocto/Shape defines functions to subdivide shape elements linearly, in
order to obtain higher shape resolution, for example before applying
displacement mapping. All functions will split all shape elements,
regardless of their size. This ensures that meshes have no cracks.
Use `subdivide_lines(lines, vert, level)` for lines,
`subdivide_triangles(triangles, vert, level)` for triangles,
`subdivide_quads(quads, vert, level)` for quads, and
`subdivide_bezier(beziers, vert, level)` for Bezier segments.
In this subdivision, each line is split in two lines,
each triangle in three triangles,
each quad in four quads, and each Bezier segment in two segments.
The functions apply the subdivision rules `level` number of times
and act on a single vertex property at a time for maximum flexibility.

```cpp
auto triangles = vector<vec3i>{...};   // initial shape
auto positions = vector<vec3f>{...};
// subdivide the triangle mesh recursively two times
auto [striangles, spositions] = subdivide_triangles(triangles, positions, 2);
```

Yocto/Shape also supports Catmull-Clark subdivision surfaces with
`subdivide_catmullclark(quads, vert, level, creased)`. In this case,
Catmull-Clark subdivision rules are used to smooth the mesh after
linear subdivision. The boundary can be treated as creases with `creased`,
which is necessary when subdividing texture coordinates.

```cpp
auto quads = vector<vec4i>{...};     // initial shape
auto positions = vector<vec3f>{...};
// subdivide the quad mesh recursively two times
auto [squads, spositions] = subdivide_catmullclark(quads, positions, 2, false);

auto tquads = vector<vec4i>{...};    // face-varying shape with texture coords
auto texcoords = vector<vec2f>{...};
// subdivide the triangle mesh recursively two times
auto [stquads, stexcoords] = subdivide_catmullclark(tquads, texcoords, 2, true);
```

## Shape sampling

Yocto/Shape supports sampling meshes uniformly. All sampling require to first
compute the shape CDF and then use it to sample the shape. For each shape type,
the sampling functions return the shape element id and the element barycentric
coordinates. Use `sample_lines(cdf, re, rn)` to sample lines,
`sample_triangles(cdf, re, rn)` to sample triangles,
`sample_quads(cdf, re, rn)` to sample quads. The shape CDFs are computed using
`sample_lines_dcf(lines, positions)`,
`sample_triangles_dcf(triangles, positions)`,
and `sample_quads_dcf(quads, positions)`.

```cpp
auto triangles = vector<vec3i>{...};   // initial shape
auto positions = vector<vec3f>{...};
auto cdf = sample_triangles_cdf(triangles, positions); // shape cdf
for(auto sample : range(samples)) {
   // sample the shape returning element id and uvs
   auto [triangle_id, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
}
```

For triangles and quads, Yocto/Shape defines convenience functions that
generate a set of points on the shape surface. Use `sample_triangles(...)`
and `sample_quads(...)` for triangles and quads respectively.

```cpp
auto triangles = vector<vec3i>{...};   // initial shape
auto positions = vector<vec3f>{...};
auto normals   = vector<vec3f>{...};
auto texcoords = vector<vec2f>{...};
auto sampled_positions = vector<vec3f>{...}; // sampled points
auto sampled_normals   = vector<vec3f>{...};
auto sampled_texcoords = vector<vec2f>{...};
// sample a set of npoints on the mesh
auto npoints = 100;
sample_triangles(sampled_positions, sampled_normals, sampled_texcoords,
                 triangles, positions, normals, texcoords, npoints);
```

## Procedural shapes

Yocto/Shape defines several procedural shapes used for both testing and
to quickly create shapes for procedural scenes. Procedural shapes take as
input the desired shape resolution, the shape scale, the uv scale, and
additional parameters specific to that procedural shape. These functions
return quads indices and vertex positions, normals and texture coordinates,
with arrays passed in.
Use `make_rect(...)` for a rectangle in the XY plane,
`make_bulged_rect(...)` for a bulged rectangle,
`make_recty(...)` for a rectangle in the XZ plane,
`make_bulged_recty(...)` for a bulged rectangle in the XZ plane,
`make_box(...)` for a box,
`make_rounded_box(...)` for a rounded box,
`make_floor(...)` for a floor in the XZ plane,
`make_bent_floor(...)` for a bent floor,
`make_sphere(...)` for a sphere obtained from a cube,
`make_uvsphere(...)` for a sphere tessellated along its uvs,
`make_capped_uvsphere(...)` for a sphere with flipped caps,
`make_disk(...)` for a disk obtained from a quad,
`make_bulged_disk(...)` for a bulged disk,
`make_uvdisk(...)` for a disk tessellated along its uvs,
`make_uvcylinder(...)` for a cylinder tessellated along its uvs,
`make_rounded_uvcylinder(...)` for a rounded cylinder.

```cpp
// most procedural shapes return quads, positions, normals, and texcoords
auto quads = vector<vec4i>{};
auto positions = vector<vec3f>{};
auto normals = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
// make shapes with 32 steps in resolution and scale of 1
make_rect(quads, positions, normals, texcoords, {32,32}, {1,1});
make_bulged_rect(quads, positions, normals, texcoords, {32,32}, {1,1});
make_recty(quads, positions, normals, texcoords, {32,32}, {1,1});
make_box(quads, positions, normals, texcoords, {32,32,32}, {1,1,1});
make_rounded_box(quads, positions, normals, texcoords, {32,32,32}, {1,1,1});
make_floor(quads, positions, normals, texcoords, {32,32}, {10,10});
make_bent_floor(quads, positions, normals, texcoords, {32,32}, {10,10});
make_sphere(quads, positions, normals, texcoords, 32, 1);
make_uvsphere(quads, positions, normals, texcoords, {32,32}, 1);
make_capped_uvsphere(quads, positions, normals, texcoords, {32,32}, 1);
make_disk(quads, positions, normals, texcoords, 32, 1);
make_bulged_disk(quads, positions, normals, texcoords, 32, 1);
make_uvdiskm(quads, positions, normals, texcoords, {32,32}, 1);
make_uvcylinder(quads, positions, normals, texcoords, {32,32,32}, {1,1});
make_rounded_uvcylinder(quads, positions, normals, texcoords, {32,32,32}, {1,1});
```

Yocto/Shape defines a few procedural face-varying shapes with similar interfaces
to the above functions. In this case, the functions return face indices and
vertex data for positions, normals and texture coordinates packed in a
`quads_fvshape` struct.
Use `make_fvrect(...)` for a rectangle in the XY plane,
`make_fvbox(...)` for a box,
`make_fvsphere(...)` for a sphere obtained from a cube.

```cpp
// procedural face-varying shapes return positions, normals, and texcoords
auto quadspos = vector<vec4i>{};
auto quadsnorm = vector<vec4i>{};
auto quadstexcoord = vector<vec4i>{};
auto positions = vector<vec3f>{};
auto normals = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
// make face-varying shapes with 32 steps in resolution and scale of 1
make_fvrect(quadspos, quadsnorm, quadstexcoord,
            positions, normals, texcoords, {32,32}, {1,1});
make_fvbox(quadspos, quadsnorm, quadstexcoord,
           positions, normals, texcoords, {32,32,32}, {1,1,1});
make_fvsphere(quadspos, quadsnorm, quadstexcoord,
              positions, normals, texcoords, 32, 1);
```

Yocto/Shape provides functions to create predefined shapes helpful in testing.
These functions take only a scale and often provide only the positions as
vertex data. These functions return either triangles, quads, or
face-varying quads.
Use `make_monkey(...)` for the Blender monkey as quads and positions only,
`make_quad(...)` for a simple quad,
`make_quady(...)` for a simple quad in the XZ plane,
`make_cube(...)` for a simple cube as quads and positions only,
`make_fvcube(...)` for a simple face-varying unit cube,
`make_geosphere(...)` for a geodesic sphere as triangles and positions only.

```cpp
// built-in shapes return elemeents, positions, normals, and texcoords
auto quads = vector<vec4i>{};
auto triangles = vector<vec3i>{};
auto quadspos = vector<vec4i>{};
auto quadsnorm = vector<vec4i>{};
auto quadstexcoord = vector<vec4i>{};
auto positions = vector<vec3f>{};
auto normals = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
// make built-in shapes with scale of 1
make_monkey(quads, positions, normals, texcoords, 1);
make_quad(quads, positions, normals, texcoords, 1);
make_quady(quads, positions, normals, texcoords, 1);
make_cube(quads, positions, normals, texcoords, 1);
make_geosphere(triangles, positions, normals, texcoords, 1);
make_fvcube(quadspos, quadsnorm, quadstexcoord,
              positions, normals, texcoords, 1);
```

Yocto/Shape supports the generation of points and lines sets.
Use `make_lines(...)` to create a line set in the XY plane,
`make_points(...)` for a collection of points at the origin,
adn `make_random_points(...)` for a point set randomly placed in a box.
These functions return shapes that are defined in terms of lines or points
and return lines or points indices, and vertex positions,
normals, texture coordinates and radia, packed in a `lines_shape`
or `points_shape` struct.

```cpp
// procedural lines return lines, positions, normals, texcoords, radia
auto lines = vector<vec2i>{};
auto positions = vector<vec3f>{};
auto normals = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto radius = vector<float>{};
make_lines(lines, positions, normals, texcoords, radius,
           {4, 65536},      // line steps and number of lines
           {1, 1}, {1, 1},  // line set scale and uvscale
           {0.001, 0.001}); // radius at the bottom and top
// procedural points return points, positions, normals, texcoords, radia
make_points(points, positions, normals, texcoords, radius,
            65536,        // number of points
            1,            // uvscale
            0.001);       // point radius
make_random_points(points, positions, normals, texcoords, radius,
            65536, // number of points
            {1, 1, 1}, 1, // line set scale and uvscale
            0.001);       // point radius
```

Yocto/Shape also defines a simple function to generate randomized hairs
on a triangle or quad mesh. Use `make_hair(...)` to create a hair shape
from a triangle and quad mesh, and return a line set.

```cpp
// Make a hair ball around a shape
auto lines = vector<vec2i>{};
auto positions = vector<vec3f>{};
auto normals = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto radius = vector<float>{};
make_hair(lines, positions, normals, texcoords, radius,
  surface_triangles, surface_quads, // sampled surface
  surface_positions, surface_normals, surface_texcoords
  {8, 65536},     // steps: line steps and number of lines
  {0.1, 0.1},     // length: minimum and maximum length
  {0.001, 0.001}, // radius: minimum and maximum radius from base to tip
  {0, 10},        // noise: noise added to hair (strength/scale)
  {0, 128},       // clump: clump added to hair (strength/number)
  {0, 0});        // rotation: rotation added to hair (angle/strength)
```

Finally, Yocto/Shape defines a function to create a quad mesh from a heighfield.
Use `make_heightfield(...)` to create a heightfield meshes.

```cpp
auto quads     = vector<vec4i>{};   // shape element buffer
auto positions = vector<vec3f>{};   // vertex data buffers
auto normals   = vector<vec3f>{};
auto texcoords = vector<vec2f>{};
auto size      = vec2i{512, 512};   // heightfield size
auto heightfield = vctor<float>{...};  // heightfield data
make_heightfield(quads, positions, normals, texcoords,
   size, heightfield);               // make heightfield mesh
```
