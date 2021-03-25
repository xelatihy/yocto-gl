# Yocto/Shading: Shading routines

Yocto/Shading defines shading and sampling functions useful to write path
tracing algorithms. Yocto/Shading is implemented in `yocto_shading.h`.

## Shading utilities

Yocto/Shading defines many functions to handle Fresnel effects, both
approximate and correct ones.
Use `fresnel_dielectric(eta, normal, outgoing)` or
`fresnel_conductor(eta, etak, normal,outgoiong)`
to evaluate the fresnel term for dielectrics and conductors,
or use `fresnel_schlick(reflectivity, normal, outgoiong)`
for the Schlick fresnel approximation.
The Schlick approximation uses reflectivity while the other functions
use the index of refraction. All these functions that the normal and
outgoing, or incoming, direction as input.

Use `eta_to_reflectivity(eta)` and `reflectivity_to_eta(reflectivity)`
to convert the dielectric index of refraction to reflectivity and vice-versa.
Use `eta_to_edgetint(eta, etak)` and `edgetint_to_eta(reflectivity, edgetint)`
to convert between the conductor index of refraction and its
artist-friendly parametrization.
Use `conductor_eta(name)` to get tabulated index of refractions for
conductors.

```cpp
// approximate Fresnel for dielectric and conductors
auto reflectivity1 = vec3f{0.04,0.04,0.04};  // dielectric reflectivity
auto reflectivity2 = vec3f{0.80,0.70,0.90};  // conductor reflectivity
auto fs1 = fresnel_schlick(reflectivity_d, normal, outgoing);
auto fs1 = fresnel_schlick(reflectivity_c, normal, outgoing);
// accurate Fresnel for dielectrics
auto ior1 = vec3f{1.5,1.5,1.5};              // dielectric ior
auto fd1 = fresnel_dielectric(ior_d, normal, outgoing);
auto ior2 = reflectivity_to_eta(reflectivity1); // reflectivity to ior
// accurate Fresnel for conductors
          {"Au", {{0.1431189557f, 0.3749570432f, 1.4424785571f},
                     {3.9831604247f, 2.3857207478f, 1.6032152899f}}},

auto eta1  = vec3f{0.143, 0.375, 1.442};     // conductor complex ior
auto etak1 = vec3f{3.983, 2.386, 1.603};     // conductor complex ior
auto fc1 = fresnel_conductor(eta1, etak1, normal, outgoing);
auto [eta2, etak2] = edgetint_to_eta(reflectivity2,vec3f{...}); // refl. to ior
auto [eta3, etak3] = conductor_eta("Au");    // get gold eta
```

Most shaders today use some version of microfacet lobes for reflections and
transmissions. Yocto/Shading has several functions to simplify writing
those shaders.
Use `microfacet_distribution(roughness, normal, halfway, ggx)` to
evaluate the microfacet distribution and
`microfacet_shadowing(roughness, normal, halfway, ggx)` to evaluate the
shadow-masking term.
Both functions take the surface roughness, the shading normal, the halfway
vector and whether or not to use the GGX or Beckman distribution.

Yocto/Shading also supports generating directions according to the microfacet
distribution, both for GGX and Beckman.
Use `sample_microfacet(roughness, normal, rn, ggx)` to generate halfway
vector with distribution proportional to the microfacet distribution
function and `sample_microfacet(roughness, normal, outgoing, rn, ggx)`
to generate dictions using the distribution of visible normals.
Use `sample_microfacet_pdf(roughness, normal, halfway, ggx)` or
`sample_microfacet_pdf(roughness, normal, halfway, outgoing, ggx)` to
compute the PDF for each sampling routine.

```cpp
auto normal   = vec3f{...};                    // shading normal
auto incoming = vec3f{...}                     // incoming direction
auto outgoing = vec3f{...};                    // outgoing direction
auto halfway = normalize(incoming, outgoing);  // halfway vector reflection
auto D = microfacet_distribution(roughness, normal, halfway); // distribution
auto G = microfacet_shadowing(roughness, normal, halfway);    // shadowing
auto halfway1 = sample_microfacet(roughness, normal, rand2f(rng)); // sample
auto pdf = sample_microfacet_pdf(roughness, normal, halfway1);// sample pdf
```

Yocto/Shading defines functions to simplify the implementation of volumetric
effects. Use `eval_transmittance(density, distance)` to evaluate the
transmittance of a homogeneous volume,
`sample_transmittance(density, max_distance, rn)` to sample a distance inside
an homogeneous medium,
`sample_transmittance_pdf(density, distance, max_distance)`
to compute the pdf of the sampled distance.
The sampling function may return a distance equal to max distance, in which
case the ray has exited the medium.
Use `mfp_to_transmission(mfp, depth)` to convert mean-free-path to transmission
at a specific depth.

```cpp
auto density = vec3f{0.99,0.99,0.99}; // medium density
auto max_distance = float{1};         // maximum distance inside the medium
auto distance = sample_transmittance(density, max_distance, // sample distance
  rand1f(rng),rand1f(rng));
auto pdf = sample_transmittance_pdf(density, distance, max_distance) // pdf
auto transmittance = eval_transmittance(density, distance);
```

Use `eval_phasefunction(anisotropy, outgoing, incoming)` to evaluate the
HG phase function,
`sample_phasefunction(anisotropy, outgoing, rn)` to sample a direction
according to the phse function,
`sample_phasefunction_pdf(anisotropy, outgoing, incoming)`
to compute the pdf of the sampled direction.

```cpp
auto outgoing = vec3f{...};           // outgoing direction
auto anisotropy = float{0};           // isotropic phase function
auto incoming = sample_phasefunction(anisotropy, outgoing, rand2f(rng));
auto phasefunc = eval_phasefunction(anisotropy, outgoing, incoming);
auto pdf = sample_phasefunction_pdf(anisotropy, outgoing, incoming);
```

## Surface materials

Yocto/Shading provides implementation for several material types.
Use `eval_<material>(<params>, normal, outgoing, incoming)` to evaluates the
lobe BSDF multiplied by the cosine,
`sample_<material>(<params>, normal, outgoing, rn)` to sample an incoming
direction, and
`sample_<material>_pdf(<params>, normal, outgoing, incoming)` to compute
the sampled direction pdf.
Yocto/Shading supports the following materials:

- `matte`: matte appearance implemented as a diffuse bsdf
- `glossy`: glossy appearance implemented as a sum of diffuse and microfacet bsdfs
- `metallic`: metallic appearance implemented as a delta or microfacet bsdfs
- `transparent`: thin glass-like appearance brdf implemented as a delta or microfacet bsdf
- `refractive`: glass-like appearance implemented as a delta or microfacet bsdf
- `passthrough`: used in volume rendering to simulated the absence of an interface
- `gltfpbr`: the pbr model used in Khronos glTF

```cpp
auto normal   = vec3f{...};            // shading normal
auto incoming = vec3f{...}             // incoming direction
auto outgoing = vec3f{...};            // outgoing direction
auto color = vec3f{1, 0.5, 0.5};       // surface color
auto roughness = float{0.1};           // roughness
auto ior = float{1.5};                 // dielectric ior
// evaluate smooth lobes
auto b1 = eval_matte(color, normal, outgoing, incoming);
auto b3 = eval_glossy(color, ior, roughness, normal, outgoing, incoming);
auto b4 = eval_metallic(color, roughness, normal, outgoing, incoming);
auto b5 = eval_transparent(color, ior, roughness, normal, outgoing,  incoming);
auto b6 = eval_refractive(color, ior, roughness, normal, outgoing, incoming);
// sample smooth lobes
auto incoming1 = sample_matte(color, normal, outgoing, rand2f(rng));
auto pdf1 = sample_matte_pdf(color, normal, outgoing, incoming)
// eval and sample delta lobes
auto incoming2 = sample_metallic(color, normal, outgoing);
auto b7 = eval_metallic(color, normal, outgoing, incoming);
```

## Design considerations

Yocto/Shading evolved from using sum of Bsdf lobes to use full Bsdfs.
Coincidentally, this evolution is similar to the way PBRT evolved from
version 3 to version 4. But the motivation for this is different.
In PBRT, the motivation for removing the sum of lobes representation was that
a sum of lobes cannot reproduce materials in a physically-correct manner.
In Yocto/Shading, the motivation is simplicity, which is a guiding principle in
Yocto/GL. In the end, all materials are just _approximations_ of
real-world behavior. So saying that a material does not properly follow
physics is weak in Yocto/GL, but entirely reasonable in systems that have
different goals. On the other hand, the resulting code is significantly simpler
and produces th same same results visually for all scenes that we have.
