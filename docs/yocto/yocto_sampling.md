# Yocto/Sampling: Sampling routines

Yocto/Sampling provides many functions to generate points and directions
useful in path tracing and procedural generation. We also include a random
number generator suitable for ray tracing.

## Random number generation

Yocto/Sampling includes an implementation of the [PCG32](https://www.pcg-random.org)
random number generator, that is a portable generator well suited for
graphics applications. The state of the generator is stored in `rng_state`
and uses only 16 bytes. The generator is default-initialized to be able to
provide random numbers as is. You can use `make_rng(seed,seq)` to initialize
a generator with a given seed and byt selecting a specific random sequence.
See [PCG32](https://www.pcg-random.org) for a discussion. If you do not need
to select a sequence, just use the default value.

Use `rand1f(rng)`, `rand2f(rng)`, `rand3f(rng)` and `rand4f(rng)` to
generate random 1-4 dimensional float vectors with coordinates in [0,1).
Use `rand1i(rng,n)` to generate a random integer in the [0,n) range.
Use `shuffle(sequence, rng)` to randomly shuffle a sequence.

```cpp
auto rng = make_rng(172784);                   // seed the generator
auto r1 = rand1f(rng); auto r2 = rand2f(rng);  // 1-2 dim. random numbers
auto r3 = rand3f(rng); auto r4 = rand4f(rng);  // 3-4 dim. random numbers
auto ri = rand1i(rng,10);                      // random int in [0,10)
auto vec = std::vector<float>{...};
shuffle(vec, rng);                             // random shuffle of a vector
```

## Generating points and directions

Yocto/Sampling defines several functions to generate random points and
directions. Each of these functions is a warp that takes random numbers in the
[0,1) domain and warps them to the desired distributions. The functions
are particularly useful in procedural modeling and path tracing. For each
function, Yocto/Sampling provides both the warp as well as its pdf, that can
be used in Monte Carlo integration.

Yocto/Sampling supports generating directions uniformly over the unit hemisphere
and sphere, or cosine distributed over the hemisphere, using respectively
`sample_hemisphere(ruv)`, `sample_sphere(ruv)`, `sample_hemisphere_cos(ruv)`
`sample_hemisphere_cos(ruv)` and `sample_hemisphere_cos(ruv)`.
For each of these functions, the corresponding `sample_<distribution>_pdf(dir)`
computes the pdf that the given direction was chosen.

```cpp
auto rng = make_rng(172784);               // seed the generator
auto rhu = sample_hemisphere(rand2f(rng)); // random hemispherical direction
auto phu = sample_hemisphere_pdf(rand2f(rng)); // direction pdf
auto rhc = sample_hemisphere_cos(rand2f(rng)); // cos-distributed direction
auto rsu = sample_sphere(rand2f(rng));     // random spherical direction
auto rhu2 = sample_hemisphere(normal,rand2f(rng));     // oriented hemisphere
auto rhc2 = sample_hemisphere_cos(normal,rand2f(rng)); // oriented hemisphere
```

Yocto/Sampling supports generating points uniformly on geometric primitives.
Use `sample_disk(uv)` to uniformly sample a disk and `sample_triangle(uv)`
to uniformly sample a triangle. For triangles we also support direct
sampling of triangle points with `sample_triangle(p0,p1,p2,ruv)`.
by warping
random numbers from the unit square to  
 uniformly over the unit hemisphere
and sphere, or cosine distributed over the hemisphere, using respectively
`sample_hemisphere(ruv)`, `sample_sphere(ruv)`, `sample_hemisphere_cos(ruv)`
`sample_hemisphere_cos(ruv)` and `sample_hemisphere_cos(ruv)`.
For each of these functions, the corresponding `sample_<distribution>_pdf(dir)`
computes the pdf that the given direction was chosen.

```cpp
auto rng = make_rng(172784);               // seed the generator
auto sdu = sample_disk(rand2f(rng));       // uniform uv on disk
auto stu = sample_triangle(rand2f(rng));   // uniform uv on triangle
auto stp = sample_triangle(p0,p1,p2,rand2f(rng));   // uniform triangle point
```

## Sampling distributions

Yocto/Sampling defines several functions to sample distributions.
Each of these functions is a warp that takes a random number in the
[0,1) domain and warps it to the desired distribution. For each
function, Yocto/Sampling provides both the warp as well as its pdf, that can
be used in Monte Carlo integration.

Use `sample_uniform(n,r)` to randomly pick an index in the `[0,n)` range,
`sample_uniform(elems,r)` to randomly pick an element from a sequence,
`sample_discrete_cdf(cdf,r)` to pick an index from a discrete distribution
that has CDF `cdf` and `sample_discrete_weights(w,r)` to pick an index from a
discrete distribution with probability equal to the weights array. For
discrete distribution, `sample_discrete_cdf()` is significantly faster when
computing many samples.

```cpp
auto rng = make_rng(172784);                 // seed the generator
auto idx = sample_uniform(num,rand1f(rng));  // uniform index
auto& elem = sample_uniform(elems,rand1f(rng));  // uniform element
auto prob = std::vector<float>{0.3,0.3,0.4};  // probability distribution
auto pidx = sample_discrete_weights(prob,rand1f(rng)); // index with prob
auto ppidx = sample_discrete_weights_pdf(prob,pidx);   // index pdf
auto cdf = std::vector<float>{prob.size()};
for(auto i : range(prob.size())               // compute cdf
  cdf[i] = prob[i] + (i ? cdf[i-1] : 0);
auto cidx = sample_discrete_cdf(cdf,rand1f(rng)); // index with cdf
auto pcidx = sample_discrete_cdf_pdf(cdf,cidx);   // index pdf
```
