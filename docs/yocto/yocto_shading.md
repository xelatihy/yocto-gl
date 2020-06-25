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

## Bsdf lobes

Most path tracing shaders are defined as sum of single BSDF lobes.
Yocto/Shading provides and implementation for several common BSDF lobes.
Three functions are implemented for each lobe. Use
`eval_<lobe>(<params>, normal, outgoing, incoming)` to evaluates the
lobe BSDF multiplied by the cosine,
`sample_<lobe>(<params>, normal, outgoing, rn)` to sample an
incoming direction, and
`sample_<lobe>_pdf(<params>, normal, outgoing, incoming)` to compute
the sampled direction pdf.
Yocto/Shading supports the following lobes:

- `diffuse_reflection`: diffuse brdf
- `diffuse_transmission`: translucency brdf, as a flipped diffuse
- `microfacet_reflection`: specular brdf for dielectrics and metals
- `microfacet_transmission`: transmission brdf for thin dielectrics
- `microfacet_refraction`: refraction brdf for dielectrics,
  includes reflection and refraction
- `delta_reflection`: delta specular for dielectrics and metals
- `delta_transmission`: delta transmission for thin dielectrics
- `delta_refraction`: delta refraction for dielectrics

```cpp
auto normal   = vec3f{...};            // shading normal
auto incoming = vec3f{...}             // incoming direction
auto outgoing = vec3f{...};            // outgoing direction
auto rs = float{0.1};                  // roughness
auto ior = float{1.5};                 // dielectric ior
auto [eta, etak] = conductor_eta("Au");// conductor complex ior
// evaluate smooth lobes
auto b1 = eval_diffuse_reflection(normal, outgoing, incoming);
auto b2 = eval_diffuse_transmission(normal, outgoing, incoming);
auto b3 = eval_microfacet_reflection(ior, rs, normal, outgoing, incoming);
auto b4 = eval_microfacet_reflection(eta,etak, rs, normal, outgoing, incoming);
auto b5 = eval_microfacet_transmission(ior, rs, normal, outgoing,  incoming);
auto b6 = eval_microfacet_refraction(ior, rs, normal, outgoing, incoming);
// sample smooth lobes
auto incoming1 = sample_diffuse_reflection(normal, outgoing, rand2f(rng));
auto pdf1 = sample_diffuse_reflection_pdf(normal, outgoing, incoming)
// eval and sample delta lobes
auto incoming2 = sample_delta_reflection(ior, normal, outgoing);
auto b7 = eval_delta_reflection(ior, normal, outgoing, incoming);
```

## Volumetric lobes

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
