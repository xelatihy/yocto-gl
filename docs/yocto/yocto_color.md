# Yocto/Color: Basic color utilities for graphics applications.

Yocto/Color provides basic color utilities for writing graphics applications,
including color color conversions, byte to float color conversions,
tone mapping and tonal adjustment.
Yocto/Color is implemented in `yocto_color.h`.

## Color Representation

Yocto/Color follows the conventions of using generic vector types
to represent color quantities rather than defining specific types,
and performing color computation in floats for increased precision.
Colors are represented as `vec3f` for three channels and `vec4f`
for four channel types. This way, colors support all arithmetic operations
defined in [Yocto/Math](yocto_math.md).

Yocto/Color supports storing colors in 8-bit representations as `vec3b`
and `vec4b`. Conversion between float and byte representation are performed
with `byte_to_float(c8)` and `float_to_byte(cf)`.

```cpp
auto red = vec3f{1,0,0}, green = vec3f{0,1,0};  // red and green colors
auto yellow = red + green, darker = red * 0.5f; // color arithmetic
auto reddish = lerp(red, green, 0.1f);          // color interpolation
auto red8 = float_to_byte(red);                 // 8bit conversion
```

## Color Conversions

## Tonal Adjustment
