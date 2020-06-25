# Yocto/Shading: Shading routines

Yocto/Shading defines shading and sampling functions useful to write path
tracing algorithms. Yocto/Shading is implemented in `yocto_shading.h`.

## Shading utilities

// 1. use `fresnel_dielectric()` or `fresnel_conductor()` to evaluate the
// fresnel term for dielectrics or conductors; use `fresnel_schlick()` for
// the Schlick fresnel approximation
// 2. use `eta_to_reflectivity()` and `reflective_to_eta()` to convert eta to
// reflectivity and vice-versa; use `eta_to_edgetint()` and
// `edgetint_to_eta()`
// 3. use `microfacet_distribution()` and `microfacet_shadowing()` to evaluate
// a microfacet distribution and its associated shadowing term

// Schlick approximation of the Fresnel term.
inline vec3f fresnel_schlick(
const vec3f& specular, const vec3f& normal, const vec3f& outgoing);
// Compute the fresnel term for dielectrics.
inline float fresnel_dielectric(
float eta, const vec3f& normal, const vec3f& outgoing);
// Compute the fresnel term for metals.
inline vec3f fresnel_conductor(const vec3f& eta, const vec3f& etak,
const vec3f& normal, const vec3f& outgoing);

// Convert eta to reflectivity
inline vec3f eta_to_reflectivity(const vec3f& eta);
// Convert reflectivity to eta.
inline vec3f reflectivity_to_eta(const vec3f& reflectivity);
// Convert conductor eta to reflectivity.
inline vec3f eta_to_reflectivity(const vec3f& eta, const vec3f& etak);
// Convert eta to edge tint parametrization.
inline std::pair<vec3f, vec3f> eta_to_edgetint(
const vec3f& eta, const vec3f& etak);
// Convert reflectivity and edge tint to eta.
inline std::pair<vec3f, vec3f> edgetint_to_eta(
const vec3f& reflectivity, const vec3f& edgetint);

// Evaluates the microfacet distribution.
inline float microfacet_distribution(float roughness, const vec3f& normal,
const vec3f& halfway, bool ggx = true);
// Evaluates the microfacet shadowing.
inline float microfacet_shadowing(float roughness, const vec3f& normal,
const vec3f& halfway, const vec3f& outgoing, const vec3f& incoming,
bool ggx = true);

// Samples a microfacet distribution.
inline vec3f sample_microfacet(
float roughness, const vec3f& normal, const vec2f& rn, bool ggx = true);
// Pdf for microfacet distribution sampling.
inline float sample_microfacet_pdf(float roughness, const vec3f& normal,
const vec3f& halfway, bool ggx = true);

// Samples a microfacet distribution with the distribution of visible normals.
inline vec3f sample_microfacet(float roughness, const vec3f& normal,
const vec3f& outgoing, const vec2f& rn, bool ggx = true);
// Pdf for microfacet distribution sampling with the distribution of visible
// normals.
inline float sample_microfacet_pdf(float roughness, const vec3f& normal,
const vec3f& halfway, const vec3f& outgoing, bool ggx = true);

## Bsdf lobes

// 4. evaluate BRDF lobes with
// - `eval_diffuse_reflection()`: diffuse brdf
// - `eval_microfacet_reflection()`: specular brdf for dielectrics and metals
// - `eval_microfacet_transmission()`: transmission brdf for thin dielectrics
// - `eval_microfacet_refraction()`: refraction brdf for dielectrics
// 5. sample BRDF lobes with `sample_XXX()` using the above lobe names
// 6. compute the PDF for BRDF lobe sampling with `sample_XXX_pdf()` using the
// above lobe names

// Evaluates a diffuse BRDF lobe.
inline vec3f eval_diffuse_reflection(
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluate a translucent BRDF lobe.
inline vec3f eval_diffuse_transmission(
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a specular BRDF lobe.
inline vec3f eval_microfacet_reflection(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a metal BRDF lobe.
inline vec3f eval_microfacet_reflection(const vec3f& eta, const vec3f& etak,
float roughness, const vec3f& normal, const vec3f& outgoing,
const vec3f& incoming);
// Evaluates a transmission BRDF lobe.
inline vec3f eval_microfacet_transmission(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a refraction BRDF lobe.
inline vec3f eval_microfacet_refraction(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);

// Sample a diffuse BRDF lobe.
inline vec3f sample_diffuse_reflection(
const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a translucency BRDF lobe.
inline vec3f sample_diffuse_transmission(
const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a specular BRDF lobe.
inline vec3f sample_microfacet_reflection(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a metal BRDF lobe.
inline vec3f sample_microfacet_reflection(const vec3f& eta, const vec3f& etak,
float roughness, const vec3f& normal, const vec3f& outgoing,
const vec2f& rn);
// Sample a transmission BRDF lobe.
inline vec3f sample_microfacet_transmission(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a refraction BRDF lobe.
inline vec3f sample_microfacet_refraction(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn);

// Pdf for diffuse BRDF lobe sampling.
inline float sample_diffuse_reflection_pdf(
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for translucency BRDF lobe sampling.
inline float sample_diffuse_transmission_pdf(
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for specular BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for metal BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(const vec3f& eta,
const vec3f& etak, float roughness, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);
// Pdf for transmission BRDF lobe sampling.
inline float sample_microfacet_transmission_pdf(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for refraction BRDF lobe sampling.
inline float sample_microfacet_refraction_pdf(float ior, float roughness,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);

// Evaluate a delta specular BRDF lobe.
inline vec3f eval_delta_reflection(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta metal BRDF lobe.
inline vec3f eval_delta_reflection(const vec3f& eta, const vec3f& etak,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta transmission BRDF lobe.
inline vec3f eval_delta_transmission(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta refraction BRDF lobe.
inline vec3f eval_delta_refraction(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);

// Sample a delta specular BRDF lobe.
inline vec3f sample_delta_reflection(
float ior, const vec3f& normal, const vec3f& outgoing);
// Sample a delta metal BRDF lobe.
inline vec3f sample_delta_reflection(const vec3f& eta, const vec3f& etak,
const vec3f& normal, const vec3f& outgoing);
// Sample a delta transmission BRDF lobe.
inline vec3f sample_delta_transmission(
float ior, const vec3f& normal, const vec3f& outgoing);
// Sample a delta refraction BRDF lobe.
inline vec3f sample_delta_refraction(
float ior, const vec3f& normal, const vec3f& outgoing, float rnl);

// Pdf for delta specular BRDF lobe sampling.
inline float sample_delta_reflection_pdf(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta metal BRDF lobe sampling.
inline float sample_delta_reflection_pdf(const vec3f& eta, const vec3f& etak,
const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta transmission BRDF lobe sampling.
inline float sample_delta_transmission_pdf(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta refraction BRDF lobe sampling.
inline float sample_delta_refraction_pdf(float ior, const vec3f& normal,
const vec3f& outgoing, const vec3f& incoming);

## Volumetric lobes

// Convert mean-free-path to transmission
inline vec3f mfp_to_transmission(const vec3f& mfp, float depth);

// Evaluate transmittance
inline vec3f eval_transmittance(const vec3f& density, float distance);
// Sample a distance proportionally to transmittance
inline float sample_transmittance(
const vec3f& density, float max_distance, float rl, float rd);
// Pdf for distance sampling
inline float sample_transmittance_pdf(
const vec3f& density, float distance, float max_distance);

// Evaluate phase function
inline float eval_phasefunction(
float anisotropy, const vec3f& outgoing, const vec3f& incoming);
// Sample phase function
inline vec3f sample_phasefunction(
float anisotropy, const vec3f& outgoing, const vec2f& rn);
// Pdf for phase function sampling
inline float sample_phasefunction_pdf(
float anisotropy, const vec3f& outgoing, const vec3f& incoming);
