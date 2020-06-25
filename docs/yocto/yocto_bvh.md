# Yocto/Bvh: Accelerated ray-intersection and point-overlap

Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
by a two-level BVH using an internal or wrapping Embree.
Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.

**This library is experimental** and will be documented appropriately when
the code reaches stability.

<!--

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
explixitly passed on evey call, while instance data uses callbacks,
since each application has its own conventions for storing those.
To make usage more convenite, we provide `bvh_XXX_data` to hold application
data and convenience wrappers for all functions.

We support working either on the whole scene or on a single shape. In the
description below yoi will see this dual API defined.

1. build the shape/scene BVH with `make_XXX_bvh()`;
2. perform ray-shape intersection with `intersect_XXX_bvh()`
3. perform point overlap queries with `overlap_XXX_bvh()`
4. refit BVH for dynamic applications with `update_XXX_bvh`

-->
