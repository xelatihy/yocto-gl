# Yocto/Geometry: Basic geometry operations for graphics applications

Yocto/Geometry defines basic geometry operations used in graphics applications,
including computation of basic geometry quantities, ray-primitive intersection,
point-primitive distance, primitive bounds, and several interpolation functions.
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
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3 = vec3f{0,1,0};
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
auto p0 = vec3f{0,0,0}, p1 = vec3f{1,0,0}, p2 = vec3f{1,1,0}, p3 = vec3f{0,1,0};
auto tp = interpolate_triangle(p0,p1,p2,{0.3,0.5});    // triangle point
auto qp = interpolate_quad(p0,p1,p2,p3,{0.3,0.5});     // quad point
auto lp = interpolate_line(p0,p1,0.3);                 // line point
auto bp = interpolate_bezier(p0,p1,p2,p3,0.3);         // bezier point
```

## Primitive bounds

// Primitive bounds.
inline bbox3f point_bounds(const vec3f& p);
inline bbox3f point_bounds(const vec3f& p, float r);
inline bbox3f line_bounds(const vec3f& p0, const vec3f& p1);
inline bbox3f line_bounds(const vec3f& p0, const vec3f& p1, float r0, float r1);
inline bbox3f triangle_bounds(
const vec3f& p0, const vec3f& p1, const vec3f& p2);
inline bbox3f quad_bounds(
const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3);

## Ray-primitive intersections

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate)
inline bool intersect_point(
const ray3f& ray, const vec3f& p, float r, vec2f& uv, float& dist);

// Intersect a ray with a line
inline bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
float r0, float r1, vec2f& uv, float& dist);

// Intersect a ray with a triangle
inline bool intersect_triangle(const ray3f& ray, const vec3f& p0,
const vec3f& p1, const vec3f& p2, vec2f& uv, float& dist);

// Intersect a ray with a quad.
inline bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
const vec3f& p2, const vec3f& p3, vec2f& uv, float& dist);

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(
const ray3f& ray, const vec3f& ray_dinv, const bbox3f& bbox);

## Point-primitive distances

// Check if a point overlaps a position pos withint a maximum distance dist_max.
inline bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p,
float r, vec2f& uv, float& dist);

// Compute the closest line uv to a give position pos.
inline float closestuv_line(const vec3f& pos, const vec3f& p0, const vec3f& p1);

// Check if a line overlaps a position pos withint a maximum distance dist_max.
inline bool overlap_line(const vec3f& pos, float dist_max, const vec3f& p0,
const vec3f& p1, float r0, float r1, vec2f& uv, float& dist);

// Compute the closest triangle uv to a give position pos.
inline vec2f closestuv_triangle(
const vec3f& pos, const vec3f& p0, const vec3f& p1, const vec3f& p2);

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
inline bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& p0,
const vec3f& p1, const vec3f& p2, float r0, float r1, float r2, vec2f& uv,
float& dist);

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
inline bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
float r2, float r3, vec2f& uv, float& dist);

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
inline bool distance_check_bbox(
const vec3f& pos, float dist_max, const bbox3f& bbox);

// Check if two bboxe overlap.
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);
