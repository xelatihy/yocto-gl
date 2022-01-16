# Yocto/Bvh: Accelerated ray-intersection and point-overlap

Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
using a two-level BVH or wrapping Intel's Embree.
Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.

## BVH representation

Yocto/Bvh provides ray-scene intersection for points, lines, triangles and
quads accelerated by a BVH data structure. Our BVH is written for
minimal code and not maximum speed, but still gives fast-enough results.

The BVH tree is stored in a `bvh_data` struct. The tree stores an array of
nodes and an array of element indices. Each node in the tree has references
to either other nodes or elements. References are represented as indices
in the nodes or elements arrays. Nodes indices refer to the nodes array,
for internal nodes, or the element arrays, for leaf nodes.
The BVH does not store shape or scene data, which is instead passed explicitly
to all calls.
We wrap this data in `shape_bvh` and `scene_bvh` to store the BVH for
[Yocto/Scene](yocto_scene.md) shapes and scenes. To support scene intersection,
we store the shapes BVHs together with the other data.

BVH nodes contain their bounds, indices to the BVH arrays of either
primitives or internal nodes, node element type,
and the split axis. Leaf and internal nodes are identical, except that
indices refer to primitives for leaf nodes or other nodes for internal nodes.

## Building BVH

Use `make_scene_bvh(scene,highquality)` or `make_shape_bvh(shape,highquaity)`
to build a bvh for a scene or shape BVH respectively. These functions takes as
input scenes and shapes from [Yocto/Scene](yocto_scene.md). By default, the BVH
is build with a fast heuristic, that can be improved slightly by setting
`highquality` to true.

```cpp
auto scene = scene_data{...};             // make a complete scene
auto bvh = make_scene_bvh(scene);         // build a BVH
```

Use `update_shape_bvh(bvh,shape)` to update a shape BVH, and
`update_scene_bvh(bvh,scene,updated_instances,updated_shapes)` to update a scene BVH,
where we indicate the indices of the instances and shapes that have been modified.
Updating works ony for change to instance frames and shapes positions.
For changes like adding or removing elements, the BVH has to be built again.

```cpp
auto scene = scene_data{...};             // make a complete scene
auto bvh = make_scene_bvh(scene);         // build a BVH
auto shapes = update_shapes(scene);       // updates some shapes
auto instances = update_instances(scene); // updates some instances
update_scene_bvh(bvh, scene, shapes, instances); // update bvh
```

## Ray intersection

Use `intersect_scene_bvh(bvh,scene,ray)` and `intersect_shape_bvh(bvh,shape,ray)`
to compute the ray-scene and ray-shape intersection respectively, and
`intersect_instance_bvh(bvh,scene,instance,ray)` to intersect a single scene instance.
These functions return a `scene_intersection` or a `shape_intersection` that
includes a `hit` flag, the `instance` index, the shape `element` index and
shape element `uv`s and the intersection `distance`.
Results values are set only if `hit` is true.
By default Yocto/Bvh computes the closet intersection, but this can be
relaxed to accept any intersection, for shadow rays, by passing an optional flag.

```cpp
auto isec = intersect_scene_bvh(bvh,scene,ray);  // ray-scene intersection
if (isec.hit) {                              // check intersection
  handle_intersection(isec.instance,         // work on intersection data
    isec.element, isec.uv, isec.distance);
}
auto isecv = intersect_scene_bvh(bvh,scene,ray,true); // check any intersection
if (!isecv.hit) {                         // check for not intersection
  handle_unoccluded(...);                 // handle unoccluded case
}
```

## Point overlap

Use `overlap_scene_bvh(bvh,scene,position,max_distance)` and
`overlap_shape_bvh(bvh,shape,position,max_distance)` to compute the scene or
shape element closest to a given point and within a maximum distance.
Use `overlap_instance_bvh(bvh,scene,instance,position,max_distance)` to test
a single scene instance. These functions return a `scene_intersection` or a
`shape_intersection`, as the intersection ones.
By default Yocto/Bvh computes the closet element, but this can be
relaxed to accept any element, by passing an optional flag.

```cpp
auto isec = overlap_scene_bvh(bvh,scene,point,dist); // closest element
if (isec.hit) {                                  // check intersection
  handle_intersection(isec.instance,             // work on intersection data
    isec.element, isec.uv, isec.distance);
}
```

## Intel's Embree Wrapper

Yocto/Bvh provide wrappers over the Intel's Embree Raytracing Kernels that
simplify their use in the context of Yocto/GL. In particular, we provide
method to build Embree Bvhs and to perform scene-element intersection, with
as API that matched the one supported by Yocto/Bvh. Embree support is enable
by setting the `YOCTO_EMBREE` compile flag and linking to Embree's libraries.

Use `make_scene_embree_bvh(scene,highquality)` or
`make_shape_embree_bvh(shape,highquaity)`,
to build a bvh for a scene or shape BVH respectively, and
`update_shape_embree_bvh(bvh,shape)` or
`update_scene_embree_bvh(bvh,scene,updated_instances,updated_shapes)`,
to update a shape or scene BVH.

```cpp
auto scene = scene_data{...};             // make a complete scene
auto bvh = make_scene_embree_bvh(scene);  // build a BVH
auto shapes = update_shapes(scene);       // updates some shapes
auto instances = update_instances(scene); // updates some instances
update_scene_embree_bvh(bvh, scene, shapes, instances); // update bvh
```

Use `intersect_scene_embree_bvh(bvh,scene,ray)` and
`intersect_shape_embree_bvh(bvh,shape,ray)`
to compute the ray-scene and ray-shape intersection respectively, and
`intersect_instance_embree_bvh(bvh,scene,instance,ray)` to intersect a single
scene instance.

```cpp
auto isec = intersect_scene_embree_bvh(bvh,scene,ray);  // intersection
if (isec.hit) {                              // check intersection
  handle_intersection(isec.instance,         // work on intersection data
    isec.element, isec.uv, isec.distance);
}
auto isecv = intersect_scene_embree_bvh(bvh,scene,ray,true); // check
if (!isecv.hit) {                         // check for not intersection
  handle_unoccluded(...);                 // handle unoccluded case
}
```
