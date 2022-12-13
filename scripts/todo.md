# ToDo V5

## Infrastructure

- Docs

## Image examples

- move to resolution and set it as default
- check sky intensity

## Shape examples

- uniform to using subdivisions for everything
- leave a grid based option, but do everything as subdivisions
- probably do everything using transformations and leave out frame
  in the generators; otherwise include frames everywhere as the only way to
  deal with sizes and positions

## Color

- simplify color space conversions by hardcoding matrices

## Random

- implement simpler Perlin noise

## IO

- consider exception in modelio
- consider exception in pbrtio
- consider directly saving ply and obj

## Math

- vecs use c array
- vecs as arrays

## Ndarrays

- decide on extents comparing to other languages
- range with sentinel
- textures as ndarrays
- potentially specialize arrays for 1d, 2d, 3d case

## Ndspan

- provide ndspans for ndarrays after testing that it works better than
  span/vector in generic code
