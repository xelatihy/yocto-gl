# Yocto/Math: Basic math types for graphics applications

Yocto/Math provides the basic math primitives used in graphics, including
small-sized vectors, matrices, frames, quaternions, rays, bounding boxes
and their transforms.

## Small vectors, matrices and frames

We provide common operations for small vectors and matrices typically used
in graphics. In particular, we support 1-4 dimensional vectors
coordinates in float and int coordinates (`vec1f`, `vec2f`, `vec3f`, `vec4f`,
`vec1i`, `vec2i`, `vec3i`, `vec4i`).

We support 2-4 dimensional matrices (`mat2f`, `mat3f`, `mat4f`) with
matrix-matrix and matrix-vector products, transposes and inverses.
Matrices are stored in column-major order and are accessed and
constructed by column. The one dimensional version is for completeness only.

To represent transformations, most of the library facilities prefer the use
coordinate frames, aka rigid transforms, represented as `frame2f` and
`frame3f`. The structure store three coordinate axes and the origin.
This is equivalent to a rigid transform written as a column-major affine
matrix. Transform operations are faster with this representation.


## Rays and bounding boxes

We represent rays in 2-3 dimensions with `ray2f`, `ray3f`.
Each ray support initialization and evaluation.

We represent bounding boxes in 2-3 dimensions with `bbox2f`, `bbox3f`.
Each bounding box support construction from points and other bounding box.
We provide operations to compute bounds for points, lines, triangles and
quads.


## Transforms

For both matrices and frames we support transform operations for points,
vectors and directions (`transform_point()`, `transform_vector()`,
`transform_direction()`). Transform matrices and frames can be
constructed from basic translation, rotation and scaling, e.g. with
`translation_mat()` or `translation_frame()` respectively, etc.
For rotation we support axis-angle and quaternions, with slerp.

