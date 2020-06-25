# Yocto/Geometry: Geometry operations

Yocto/Geometry defines basic geometry operations, including computation of
basic geometry quantities, ray-primitive intersection, point-primitive
distance, primitive bounds, and several interpolation functions.
Yocto/Geometry is implemented in `yocto_geometry.h`.

## Primitive parametrization

Yocto/Geometry supports operations on lines, triangles, quads and beziers.
In these functions, triangles are parameterized with barycentric coordinates
`uv` written with respect to the `p1-p0` and `p2-p0` axes respectively.
Quads are internally handled as pairs of triangles `p0,p1,p3` and `p2,p3,p1`,
with the `uv` coordinates of the second triangle corrected as `1-u` and `1-v`
to produce a quad parametrization where `u` and `v` go from 0 to 1.
Degenerate quads with `p2==p3` represent triangles correctly.
This parametrization is equivalent to Intel's Embree.

## Geometric properties

For lines, Yocto/Geometry supports the computation of line lengths,
with `line_length(p0,p1)`, and line tangents, with `line_tangent(p0,p1)`.
For triangles, Yocto/Geometry supports the computation of triangle areas,
with `triangle_area(p0,p1,p2)`, and triangle normals,
`triangle_normal(p0,p1,p2)`. Similarly for quads, use `quad_area(p0,p1,p2,p3)`,
`quad_normal(p0,p1,p2,p3)` for areas and normals.

For triangles and quads, Yocto/Geometry also supports computing tangents and
bitangents from texture coordinates. This is helpful for applying normal or
bump mapping during rendering.
Use `triangle_tangents_fromuv(p0,p1,p2,uv0,uv1,uv2)` for triangles
and `quad_tangents_fromuv(p0,p1,p2,p3,uv0,uv1,uv2,uv3,uv)` for quads.

For triangles, tangents and bitangents are defined with respect to the first
vertex as the origin.
For quads, we define the vectors on the two triangles and do not compute the
average. For this pass an additional texture coordinate since internally we split
the triangle into two and we need to known where to do it.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto ta = triangle_area(p0,p1,p2);   // triangle area
auto qa = quad_area(p0,p1,p2,p3);    // quad area
auto tn = triangle_normal(p0,p1,p2); // triangle normal
auto qn = quad_normal(p0,p1,p2,p3);  // quad normal
auto uv0 = vec2f{0,0}, uv1 = vec2f{1,0}, uv2 = vec2f{1,1}, uv3 = vec2f{0,1};
auto [tu, tv] = triangle_tangents_fromuv(p0,p1,p2,uv0,uv1,uv2); // tangents
```

## Interpolation on primitives

For all primitives, Yocto/Geometry defines interpolation functions that take
values defined at the primitive vertices and compute the interpolate value
at the parametrized point. Lines and beziers are parametrized with their natural
parameter, while triangles and quads use barycentric interpolation as defined
above.

Use `interpolate_line(p0,p1,u)` for lines, `interpolate_triangle(p0,p1,p2,uv)`
for triangles, `interpolate_quad(p0,p1,p2,p3,uv)` for quads and
`interpolate_bezier(p0,p1,p2,p3,u)` for cubic Bezier segments,
whose derivatives can be computed with
`interpolate_bezier_derivative(p0,p1,p2,p3,u)`.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto tp = interpolate_triangle(p0,p1,p2,{0.3,0.5});  // triangle point
auto qp = interpolate_quad(p0,p1,p2,p3,{0.3,0.5});   // quad point
auto lp = interpolate_line(p0,p1,0.3);               // line point
auto bp = interpolate_bezier(p0,p1,p2,p3,0.3);       // bezier point
```

## Primitive bounding boxes

Yocto/Geometry provides functions to compute bounding boxes for all primitives
types. For points and lines, vertices might have a thickness associate with them.
Use `point_bounds(p0,r0)` for points, `line_bounds(p0,p1,r0,r1)` for lines,
`triangle_bounds(p0,p1,p2)` for triangles, `quad_bounds(p0,p1,p2,p3)` for quads.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto tb = triangle_bounds(p0,p1,p2);  // triangle bounding box
auto qb = quad_bounds(p0,p1,p2,p3);   // quad bounding box
auto r0 = 0.01, r1 = 0.01;
auto lb = line_bounds(p0,p1,r0,r1);   // line bounding box
auto pb = point_bounds(p0,r0);        // point bounding box
```

## Ray-primitive intersections

Yocto/Geometry defines functions for ray-primitive intersection. Each function
returns wether the primitive was hit and, if so, sets the primitive parameters
and the intersection distance as output variables. Triangle intersection are
computed using the Moller-Trombone intersection algorithm.
Quad intersections are computed
by treating quads as two triangles. Point intersections are compute
approximately, by treating points as ray-oriented disks. Line intersections
are computed approximately, by treating lines as ray-oriented ribbons.
Use `intersect_point(ray,p0,r0,uv,d)` for points,
`intersect_line(ray,p0,p1,r0,r1,uv,d)` for lines,
`intersect_triangle(ray,p0,p1,p2,uv,d)` for triangles,
`intersect_quad(ray,p0,p1,p2,p3,uv,d)` for quads.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto r0 = 0.01, r1 = 0.01;
auto ray = ray3f{{0,0,0.5},{0,0,-1}};                 // ray
auto uv = vec2f{0,0}; auto dist = float{0};           // hit distance and uvs
auto th = intersect_triangle(ray,p0,p1,p2,uv,dist));  // triangle intersection
auto qh = intersect_quad(ray,p0,p1,p2,p3,uv,dist));   // quad intersection
auto lh = intersect_line(ray,p0,p1,r0,r1,uv,dist));   // line intersection
auto ph = intersect_point(ray,p0,r0,uv,dist));        // point intersection
```

Yocto/Geometry defines two functions to test whether a ray hits a bounding box.
In this case, we do not return the ray distance or hit, but just check for
intersection, which is useful when defining BVH hierarchies.
Use `intersect_bbox(ray,bbox)` as a simple alternative and
`intersect_bbox(ray,ray_dinv,bbox)` for a faster one.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto bbox = quad_bounds(p0,p1,p2,p3);
auto bh = intersect_bbox(ray, bbox);           // bbox intersection check
auto ray_dinv = 1 / ray.d;                     // ray direction inverse
auto bf = intersect_bbox(ray, ray_dinv, bbox); // fast bbox intersection
```

## Point-primitive overlaps

Yocto/Geometry defines functions for point-primitive distance and overlap.
Each function returns wether the primitive was hit and, if so, sets
the primitive parameters and the overlap distance as output variables.
Each function takes a position and a maximum distance to test within,
together with primitive vertices and thickness.
Use `overlap_point(pt,md,p0,r0,uv,d)` for points,
`overlap_line(pt,md,p0,p1,r0,r1,uv,d)` for lines,
`overlap_triangle(pt,md,p0,p1,p2,r0,r1,r2,uv,d)` for triangles,
`overlap_quad(pt,md,p0,p1,p2,p3,r0,r1,r2,r3,uv,d)` for quads.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto r0 = 0.01, r1 = 0.01, r2 = 0.01, r3 = 0.01;
auto pt = vec3f{0,0,0.5}; auto md = float{1};         // point and max dist
auto uv = vec2f{0,0}; auto dist = float{0};           // hit distance and uv
auto th = overlap_triangle(pt,md,p0,p1,p2,uv,dist));  // triangle overlap
auto qh = overlap_quad(pt,md,p0,p1,p2,p3,uv,dist));   // quad overlap
auto lh = overlap_line(pt,md,p0,p1,r0,r1,uv,dist));   // line overlap
auto ph = overlap_point(pt,md,p0,r0,uv,dist));        // point overlap
```

Yocto/Geometry defines a function to test whether a point is contained within a
bounding bbox within a certain distance. Just like before, we do not return
the ray distance or hit, but just check for overlap, which is useful when
defining BVH hierarchies.
Use `overlap_bbox(pt,md,bbox)` to test for overlap between a point and
a bounding box and `overlap_bbox(bbox1,bbox2)` to test whether two bounding
boxes overlap.

```cpp
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3=vec3f{0,1,0};
auto bbox = quad_bounds(p0,p1,p2,p3);
auto bbox2 = bbox3f{{0,0,0}, {1,1,1}};
auto pt = vec3f{0,0,0.5}; auto md = float{1};   // point and max dist
auto bh  = overlap_bbox(pt, md, bbox);          // bbox overlap check
auto bh2 = overlap_bbox(bbox, bbox2);           // bbox overlap check
```
