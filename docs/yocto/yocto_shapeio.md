# Yocto/ShapeIO: Shape serialization

Yocto/ShapeIO supports loading and saving shapes from Ply, Obj, Stl.
Yocto/ShapeIO is implemented in `yocto_shapeio.h` and `yocto_shapeio.cpp`.

## Shape serialization

Yocto/ShapeIO supports reading and writing shapes. 
Use `load_shape(filename, shape)` or `shape = load_shape(filename)` 
to load shapes, and `save_shape(filename, shape)` to save it.
Shapes are stored as `shape_data` structs defined in [Yocto/Shape](yocto_shape.md).
Upon errors, an `io_error` is thrown from all IO functions.
See [Yocto/CommonIO](yocto_commonio.md) for discussion on error handling 
and use without exceptions.

```cpp
auto shape = shape_data{};
load_shape("input_file.ply",  shape);        // load shape
save_shape("output_file.ply", shape);        // save shape
auto shape1 = load_shape("input_file.ply");  // alternative load
```

Yocto/ShapeIO supports reading and writing face-varying shapes. 
Use `load_fvshape(filename, shape)` or `shape = load_fvshape(filename)` 
to load shapes, and `save_shape(filename, shape)` to save it.
Shapes are stored as `fvshape_data` structs defined in [Yocto/Shape](yocto_shape.md).
Upon errors, an `io_error` is thrown from all IO functions.
See [Yocto/CommonIO](yocto_commonio.md) for discussion on error handling 
and use without exceptions.

```cpp
auto shape = fvshape_data{};
load_fvshape("input_file.obj",  shape);        // load shape
save_fvshape("output_file.obj", shape);        // save shape
auto shape1 = load_fvshape("input_file.obj");  // alternative load
```
