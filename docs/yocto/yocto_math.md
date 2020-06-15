# Yocto/Math: Basic math types for graphics applications

Yocto/Math defines the basic math primitives used in graphics, including
small-sized vectors, matrices, frames, quaternions, rays, bounding boxes
and their transforms. Yocto/Math is implemented in `yocto_math.h`.

## Vectors

The most basic types defined in Yocto/Math are fixed-length vectors in
1 to 4 dimensions with either `float` or `int` coordinates, namely `vec1f`,
`vec2f`, `vec3f`, `vec4f`, `vec1i`, `vec2i`, `vec3i`, `vec4i`.
One-dimensional vectors are defined for completeness.

Vectors coordinates are default-initialized to zeros. Vectors are constructed
by specifying all their coordinates. For convenience, Yocto/Math defines
common constants, i.e. `zeroXX`, to initialize zero vectors.
Vector coordinates are accessed by index, e.g. `v[0]`, or by name, using
`v.x`, `v.y`, `v.z`, `v.w`. While in graphics it is common to use the latter,
we suggest to instead use the former since it better generalizes to
different cases.

Vectors support component-wise arithmetic, equality and assignment operators.
For arithmetic operators, both vector and scalar operands are supported.
Yocto/Math also overrides many common math function to work on vectors
by component-wise evaluation, including trigonometric and exponentiation
functions, and component-wise abs, max, min and clamp. Real vector types
additionally support vector products and lengths.

## Matrices

Yocto/Math supports 2-4 dimensional square matrices, namely `mat2f`, `mat3f`
and `mat4f`. Matrices are stored on column-major order. Matrices are
initialized to the identity matrix and constructed by passing matrix columns.
Matrix columns are accessed using indexing, i.e. `mat[col]`.
Yocto/Math defines convenience constants for common matrices, e.g.
`indentity3x3f`, etc.

Matrices support linear algebra summation and scaling operations,
component-wise equality, matrix-matrix and matrix-vector products.
All these operators use standard math notation.
Matrices also support the computation of transpose, inverse and adjoint
matrices, with `transpose(a)`, `inverse(a)` and `adjoint(a)`.

## Frames

To represent transformations, most of the library facilities prefer the use
coordinate frames, aka rigid transforms, represented as `frame2f` and
`frame3f`. Frames are represented in column-major order, with columns
defined as the axis of the coordinate frame followed by the origin.
This is equivalent to representing rigid transformations as column-major affine
matrices.

Frames support only multiplication operations with other frames, to chain
transforms, the computation of their inverses, with `inverse(f)`, and the
creation of frames aligned with one ot two axis, i.e. `frame_fromz(o,z)`
and `frame_fromxz(o,z,x)`. For the inverse computation, the default is to
assume that the frame is orthonormal to get a fast answer. This can be
overridden as shown in the example.

```cpp
auto f1 = frame3f{x, y, z, o}; // construct a frame from axes and origin
auto of = inverse(f1);         // assume f is orthonormal
auto gf = inverse(f1,true);    // treat f as a generic affine matrix
```

## Rays

Yocto/Math defines rays in 2-3 dimensions as `ray2f` and `ray3f`.
Rays are defined as an origin `o`, a direction `d` and minimum and maximum
values for the distance along a ray, namely `tmin` and `tmax`.
Besides construction, Yocto/Math does not specify other ray functionality
which is instead use in other parts of Yocto/GL.

## Bounding boxes

Yocto/Math defines axies-aligned bounding boxes in 2 to 3 dimensions as
`bbox2f` and `bbox3f`. Bounding boxes store the minimum and maximum coordinate
values, that can be accessed with `b.min` and `b.max`. Bounding boxes are
default-initialized to an invalid state that contains no points,
or they are constructed by specifying the min and max values directly.

To build bounds for complex primitives, bounding boxes are very initialized to
empty bounds, that can be done by using the constants like `invalidabXf`,
and then grown to encompass either points or other bounding boxes with
`merge(b,p)`.

```cpp
auto bbox = invalidb3f;
for(auto point : points) bbox = merge(bbox, point);
```

## Transforms

Yocto/Math uses linear algebra vectors for holding different geometric
quantities, such as points, vectors, directions and normals. For this reason,
we define different functions to transform these quantities, namely

- `transform_point(f,p)` to transform points,
- `transform_vector(f,v)` to transform vectors,
- `transform_direction(f,d)` to transform directions intended as
  normalized vectors,
- `transform_normal(f,n)` to transform normals.

Yocto/Math provides overrides to transforms for both frames and matrices,
but _the preferred transform representation in Yocto/GL are frames_.
When transforming normals, the inverse transpose is used for matrices,
while a
while frames use a fast inverse or an affine inverse as
and just
