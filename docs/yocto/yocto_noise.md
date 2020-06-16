# Yocto/Noise: Noise functions

Yocto/Noise provides a Perlin noise implementation.
_This library should to be considered a placeholder since it will grow later
as a collection of noise functions_ used in procedural modeling.
For now, the implementation used is the one
found in the [stb libraries](https://github.com/nothings/stb),
that are released in the public domain.
Yocto/Noise is implemented in `yocto_noise.h`.

## Noise functions

Use `perlin_noise(p,w)` to generate Perlin noise with optional wrapping.
For fractal variations, use `perlin_ridge(p,l,g,o,f,w)`,
`perlin_fbm(p,l,g,o,w)` and `perlin_turbulence(p,l,g,o,w)`.
Each fractal version is defined by its lacunarity `l`, its gain `g`, the
number of octaves `o` and possibly an offet.

```cpp
auto p = vec3f{0,0,0};
auto n = perlin_noise(p);
auto lacunarity = 2.0f, gain = 0.5.0f, , offset = 1.0f; auto octaves = 6;
auto n = perlin_ridge(p, lacunarity, gain, octaves, offset);
auto n = perlin_fbm(p, lacunarity, gain, octaves);
auto n = perlin_turbulence(p, lacunarity, gain, octaves);
```
