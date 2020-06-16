# Yocto/Sampling: Sampling routines for graphics applications

Yocto/Sampling provides many functions to generate points and directions
useful in path tracing and procedural generation. We also include a random
number generator suitable for ray tracing.

## Random Number Generation

This library supports generting random numbers using the PCG32 algorithm,
that is a portable generator well suited for graphics applications.

1. initialize the random number generator with `make_rng()`
2. if necessary, you can reseed the rng with `seed_rng()`
3. generate random integers in an interval with `rand1i()`
4. generate random floats and double in the [0,1) range with `rand1f()`,
   `rand2f()`, `rand3f()`, `rand1d()`

## Generating points and directions

1. use `sample_XXX()` to warp random numbers in [0,1)^k domains to the
   desired domain; in particular we support `sample_hemisphere()`,
   `sample_sphere()`, `sample_hemisphere_cos()`,
   `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
   `sample_triangle()`, `sample_quad()`
2. use `sample_discrete()` to sample from a descreet distribution
3. use `sample_XXX_pdf()` to compute the PDF of the sampling functions
