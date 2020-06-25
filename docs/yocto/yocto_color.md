# Yocto/Color: Color operations

Yocto/Color defines basic color manipulations, including color color
conversions, byte to float color conversions, tone mapping and tonal
adjustments. Yocto/Color is implemented in `yocto_color.h`.

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

Like most 3D graphics libraries, Yocto/Color does not explicitly track color spaces. Instead, we use conventions between library facilities to carry out
color operations. Most color operations are defined on linear RGB colors.
By default, Yocto/Color uses a linear color space with sRGB primaries, mostly
since lots of freely available graphics data encodes colors in that manner.

Yocto/Color supports conversions between linear RGB, sRGB, XYZ, xyY and HSV.
Color conversion functions are named as `<from>_to_<to>` where `<from>` and
`<to>` are the name of the color spaces.

```cpp
auto c_rgb = vec3f{1,0.5,0};      // color in linear RGB
auto c_xyz = rgb_to_xyz(c_rgb);   // convert to XYZ
auto c_srgb = rgb_to_srgb(c_rgb); // convert to sRGB
auto c_hsv = rgb_ro_hsv(c_srgb);  // sRGB to HSV
auto c_8bit = float_to_byte(rgb_to_srgb(c_rgb)); // encode in 8-bit sRGB
```

## Tonal Adjustment

Yocto/Color provides minimal tonal adjustment functionality, mostly helpful
when displaying images on screens. HDR images can be tone mapped with
`tonemap(hdr,e,f)` that applies a simple exposure correction to the HDR colors
and optionally a filmic tone curve. This functions converts from linear
RGB to sRGB by default, but the latter conversion can be disabled.

Saturation adjustment are implemented in `sature(rgb,s)` that interpolates
to color with gray.

Contrast adjustments are notoriously more complex since they depend on
both the contrast function used and the color space it is applied with. To
cater to most uses, Yocto/Color provides three contrast curves.
`contrast(rgb,c)` generate smooth contrast variations by applying an S-shaped
curve that typically used in non-linear color spaces.
`lincontrast(rgb,c,g)` has a harsher look due to clipped colors and is obtained
by applying linear scaling around a gray level.
`logcontrast(rgb,c,g)` is suitable to contrast manipulation in HDR since it
does not clip colors. This is implemented by applying scaling in the
logarithm domain. Both linear and logarithmic contrast can be used for HDR and
LDR by adjusting the gray level to 0.18 and 0.5 respectively.

```cpp
auto exposure = 1; auto filomic = true;
auto ldr = tonemap(hdr, exposure, filmic);        // HDR to LDR tone mapping
ldr = contrast(saturate(ldr,0.7),0.7);            // contrast and saturation
auto hdr2 = logcontrast(hdr * tint, 0.7, 0.18);   // contrast and tint in HDR
```
