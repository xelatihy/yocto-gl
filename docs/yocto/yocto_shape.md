# Yocto/Shape: Shape utilities

**WORK IN PROGRESS**

Yocto/Shape is a collection of utilities for manipulating shapes in 3D
graphics, with a focus on triangle and quad meshes.

## Shape functions

We provide a small number of utilities for shape manipulation for index
triangle and quad meshes, indexed line and point sets and indexed beziers.
The utilities collected here are written to support a global illumination
rendering and not for generic geometry processing. We support operation for
shape smoothing, shape subdivision (including Catmull-Clark subdivs), and
example shape creation.

1. compute smooth normals and tangents with `compute_normals()`
   `compute_tangents()`
2. compute tangent frames from texture coordinates with
   `compute_tangent_spaces()`
3. compute skinning with `compute_skinning()` and
   `compute_matrix_skinning()`
4. create shapes with `make_proc_image()`, `make_hair()`,
   `make_points()`, `make_point()`
5. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
6. shape sampling with `sample_points()`, `sample_lines()`,
   `sample_triangles()`; initialize the sampling CDFs with
   `sample_points_cdf()`, `sample_lines_cdf()`,
   `sample_triangles_cdf()`
7. sample a could of point over a surface with `sample_triangles()`
8. get edges and boundaries with `get_edges()`
9. convert quads to triangles with `quads_to_triangles()`
10. convert face varying to vertex shared representations with
    `convert_face_varying()`
11. subdivide elements by edge splits with `subdivide_lines()`,
    `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
12. Catmull-Clark subdivision surface with `subdivide_catmullclark()`

## Shape IO

We support reading and writing shapes in OBJ and PLY.

1. load/save shapes with `load_shape()`/`save_shape()`

## Ray-Scene and Closest-Point Queries

Yocto/BVH is a simple implementation of ray intersection and
closest queries using a two-level BVH data structure. We also include
low-level intersection and closet point primitives.
Alternatively the library also support wrapping Intel's Embree.

Yocto/GL provides ray-scene intersection for points, lines, triangles and
quads accelerated by a two-level BVH data structure. Our BVH is written for
minimal code and not maximum speed, but still gives reasonable results. We
suggest the use of Intel's Embree as a more efficient alternative.

In Yocto/Bvh, shapes are described as collections of indexed primitives
(points/lines/triangles/quads) like the standard triangle mesh used in
real-time graphics. A scene if represented as transformed instances of
shapes. The internal data structure is a two-level BVH, with a BVH for each
shape and one top-level BVH for the whole scene. This design support
instancing for large scenes and easy BVH refitting for interactive
applications.

In these functions triangles are parameterized with uv written
w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
as pairs of two triangles p0,p1,p3 and p2,p3,p1, with the u/v coordinates
of the second triangle corrected as 1-u and 1-v to produce a quad
parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
represent triangles correctly, an this convention is used throught the
library. This is equivalent to Intel's Embree.

Shape and scene data is not copied from the application to the BVH to
improve memory footprint at the price of convenience. Shape data is
explicitly passed on every call, while instance data uses callbacks,
since each application has its own conventions for storing those.
To make usage more convenient, we provide `bvh_XXX_data` to hold application
data and convenience wrappers for all functions.

We support working either on the whole scene or on a single shape. In the
description below yoi will see this dual API defined.

1. build the shape/scene BVH with `make_XXX_bvh()`;
2. perform ray-shape intersection with `intersect_XXX_bvh()`
3. perform point overlap queries with `overlap_XXX_bvh()`
4. refit BVH for dynamic applications with `update_XXX_bvh`
