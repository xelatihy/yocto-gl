# Yocto/Math: Math types

Yocto/Math defines the basic math primitives used in graphics, including
small-sized vectors, matrices, frames and quaternions, and their transforms.
Yocto/Math is implemented in `yocto_math.h`.

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
additionally support vector products and lengths. For 3D vectors, we also
define reflection and refraction vectors following GLSL.

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
This is equivalent to representing rigid transformations as column-major
affine matrices. Use `mat_to_frame(mat)` and `frame_to_mat(frame)` to convert
between these types.

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

## Quaternions

Quaternions are represented as `quat4f`. Quaternions support is minimal and
focuses on representing rotations, i.e. unit quaternions. Quaternions are
initialized to a unit quaternion, that can also be set with the convenience
constant `identityq4f`.

Quaternions support addition, multiplication and scaling, together with
`dot(q1,q2)`, `length(q)` and `normalize(q)`. For unit quaternions,
Yocto/Math also defines the quaternion inverse `inverse(q)` and several
interpolation function including linear interpolation, normalized linear
interpolation and spherical interpolation, respectively with `lerp(a,bt)`,
`nlerp(a,b,t)` and `slerp(a,b,t)`.

## Transforms

Yocto/Math uses linear algebra vectors for holding different geometric
quantities, such as points, vectors, directions and normals. For this reason,
we define different functions to transform these quantities, namely
`transform_point(f,p)` to transform points, `transform_vector(f,v)` to
transform vectors, `transform_direction(f,d)` to transform directions
intended as normalized vectors, and `transform_normal(f,n)` to transform
normals.

Yocto/Math provides overrides to transforms for both frames and matrices,
but _the preferred transform representation in Yocto/GL are frames_.
When transforming normals by matrices, the inverse transform is computed on
the fly for correctness. When transforming normals by frames, frames are
assumed to be orthonormal and used as is, unless specifically requested as
for inverses.

For convenience, Yocto/Math provides several functions to construct transform
frames and matrices. Translation, rotation and scaling frames are created with  
`translation_frame(t)`, `rotation_frame(r)` and `scaling_frame(s)`. Rotation
frames can be derived from axis-angle, quaternion and matrix representations.
To define camera frames, one can use `lookat_frame(from,to,up)`.

```cpp
auto f = tranlation_frame({1,0,0}) * rotation_frame({0,1,0},pi/2); // transform
auto lp = vec3f{1,2,3}, lv = vec3f{0,2,3}; // point and vector in local coords
auto wp = transform_point(f, lp);          // point in world coords
auto wv = transform_vector(f, lv);         // vector in world coords
```

For use in GPU programming, Yocto/Math also defines various projection
matrices in the style of GLU/GLM, namely `frustum_mat(...)`,
`ortho_mat(...)`, `ortho2d_mat(...)`, and `perspective_mat(...)`.

## User-Interface Transforms

Yocto/Math provides a few utilities for writing user interfaces for 2D images
and 3D environments. For images, we assume that image are scaled and than
translated to be place oin screen. User interfaces can directly manipulate
the translation and zoom without further helpers. Yocto/Math just provides
convenience functions for centering an image and compute a mouse to image
coordinate transformation.

```cpp
auto image_size = vec2i{512,512}, window_size = vec2i{1024,768};
auto image_center = vec2f{100,100}; auto image_scale = float{1};
// get image coordinates from mouse
auto ij = get_image_coords(mouse_pos, image_center, image_scale, image_size);
// center image by setting image_center and image_scale
update_imview(image_center, image_scale, image_size, window_size, false);
// center image and for to screen by setting image_center and image_scale
update_imview(image_center, image_scale, image_size, window_size, true);
```

For 3D cameras, Yocto/Math supports a turntable interface, inspired by many 3D
editors, with `update_turntable(...)`, and the camera's rays generation with
`camera_ray(...)`.

```cpp
auto camera_frame = identity3x4f; auto camera_focus = float{10};
// update camera position and focus
update_turntable(camera_frame, camera_focus, rotate, dolly, pan);
// turntable udpates works also for camera with explicit lookat parametrizations
update_turntable(camera_from, camera_to, camera_up, rotate, dolly, pan);
// generation of camera rays
auto ray = camera_ray(camera_frame, lens, aspect, film, image_uv);
```
