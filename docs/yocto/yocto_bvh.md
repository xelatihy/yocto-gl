# Yocto/Bvh: Accelerated ray-intersection and point-overlap

Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
using a two-level BVH or wrapping Intel's Embree.
Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.

## BVH representation

Yocto/Bvh provides ray-scene intersection for points, lines, triangles and
quads accelerated by a BVH data structure. Our BVH is written for
minimal code and not maximum speed, but still gives fast-enough results.

The BVH tree is stored in a `bvh_tree` struct. The tree stores an array of
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

Two wrappers, `bvh_scene` and `bvh_shape` are used to store the BVH for
[Yocto/Scene](yocto_scene.md) scenes and shapes. For shapes, we store a
single BVH, while for scene we store the scene BVH and all shapes BVHs.
These wrappers does not store copies to shape or scene data, so that data is
passed in for all subsequent calls.

These wrappers can alternatively store references to Intel's Embree BVH,
for faster ray-scene intersection. To use Embree, the library should be
compiled with Embree support by setting the `YOCTO_EMBREE` compile flag and
linking to Embree's libraries.

## Building BVH

Use `make_bvh(scene,highquality,embree)` or `make_bvh(shape,highquaity,embree)`
to build a bvh for a scene or shape BVH respectively.
These functions takes as input scenes and shapes from
[Yocto/Scene](yocto_scene.md). By default, the BVH is build with a fast heuristic,
that can be improved slightly by setting `highquality` to true.
By default, Yocto/BVH uses the internal BVH. Intel's Embree can be used
by setting the `embree` flag to true.

```cpp
auto scene = scene_model{...};            // make a complete scene
auto bvh = build_bvh(scene);              // build a BVH
auto embree = build_bvh(scene,true,true); // use Embree
```

Use `update_bvh(bvh,shape)` to update a shape BVH, and
`update_bvh(bvh,scene,updated_instances,updated_shapes)` to update a scene BVH,
where we indicate the indices of the instances and shapes that have beed modified.
Updating works ony for change to instance frames and shapes positions.
For changes like adding or removing elements, the BVH has to be built again.

```cpp
auto scene = scene_model{...};            // make a complete scene
auto bvh = build_bvh(scene);              // build a BVH
auto shapes = update_shapes(scene);       // updates some shapes
auto instances = update_instances(scene); // updates some instances
update_bvh(bvh, scene, shapes, instances);// update bvh
```

## Ray intersection

Use `intersect_bvh(bvh,scene,ray)` and `intersect_bvh(bvh,shape,ray)` to
compute the ray-scene and ray-shape intersection respectively, and
`intersect_bvh(bvh,scene,instance,ray)` to intersect a single scene instance.
These functions return a `bvh_intersection` that includes a `hit` flag,
the `instance` index, the shape `element` index and shape element `uv`s and
the intersection `distance`. Results values are set only if `hit` is true.
By default Yocto/Bvh computes the closet intersection, but this can be
relaxed to accept any intersection, for shadow rays, by passing an optional flag.

```cpp
auto isec = intersect_bvh(bvh,scene,ray); // ray-scene intersection
if (isec.hit) {                           // check intersection
  handle_intersection(isec.instance,      // work on intersection data
    isec.element, isec.uv, isec.distance);
}
auto isecv = intersect_bvh(bvh,scene,ray,true); // check any intersection
if (!isecv.hit) {                         // check for not intersection
  handle_unoccluded(...);                 // handle unoccluded case
}
```

## Point overlap

Use `overlap_bvh(bvh,scene,position,max_distance)` and
`overlap_bvh(bvh,shape,position,max_distance)` to compute the scene or
shape element closest to a given point and within a maximum distance.
Use `overlap_bvh(bvh,scene,instance,position,max_distance)` to test
a single scene instance. These functions return a `bvh_intersection`
as the intersection ones.
By default Yocto/Bvh computes the closet element, but this can be
relaxed to accept any element, by passing an optional flag.

```cpp
auto isec = overlap_bvh(bvh,scene,point,dist); // closest element
if (isec.hit) {                                // check intersection
  handle_intersection(isec.instance,      // work on intersection data
    isec.element, isec.uv, isec.distance);
}
```
