# Yocto/Geometry: Basic geometry support for graphics applications.

Yocto/Geometry provides basic geometry operations used in graphics applications.

## Geometry functions

The library supports basic geometry functions such as computing
line/triangle/quad normals and areas, picking points on triangles
and the like. In these functions, triangles are parameterized with uv written
w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
as pairs of two triangles (p0,p1,p3) and (p2,p3,p1), with the uv coordinates
of the second triangle corrected as 1-u and 1-v to produce a quad
parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
represent triangles correctly, and this convention is used throught the
library. This is equivalent to Intel's Embree.
