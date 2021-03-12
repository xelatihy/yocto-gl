# Yocto/Bvh: Accelerated ray-intersection and point-overlap

Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
using a two-level BVH or wrapping Intel's Embree.
Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.

## Building BVH

Use `make_bvh(scene,highquality,embree)` or `make_bvh(shape,highquaity,embree)`
to build a bvh for a scene or shape BVH respectively.
These functions takes as input scenes and shapes from
[Yocto/Scene](yocto_scene.md). By default, the BVH is build with a fast heuristic,
that can be improved slightly by setting `highquality` to true. By default, we
use the internal BVH, but Intel's Embree can be used, when available, if
`embree` is set to true.
Yocto/BVH works in shared memory so the data structures returned do not
copy scene or shape data.

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
