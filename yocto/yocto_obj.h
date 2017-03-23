///
/// YOCTO_OBJ: Wavefront OBJ/MTL loader and writer with support for points,
/// lines, triangles and general polygons and all materials properties.
/// Contains also a few extension to eqasily create demos such as per-vertex
/// color and radius, cameras and envmaps. Can use either a low-level OBJ
/// representation or a high level flattened representation.
///
///
/// USAGE FOR READING:
///
/// 1. load an obj with load_obj()
///   - loads an obj from disk including its associate mtl files
///   - returns a parsed scene data structure described below
/// 2. [LOW-LEVEL INTERFACE] access the data directly from the returned object
///   - the data is documented below and matches the OBJ file structure exactly
/// 3. [HIGH-LEVEL INTERFACE] optionally flatten the data as a more friendly
///    representation where shapes are index meshes, supporting points, lines
///    and triangle primitives, with flatten_obj()
///   - the flattened data, documented below, can be use to draw directly on
///   the GPU or in a raytracer
///   - vertices are duplicated as needed to support GPU friendly access
///   - optionally load textures data as float arrays with
///       load_fl_textures()
///
/// Both in reading and writing, OBJ has no clear convention on the orientation
/// of textures Y axis. So in many cases textures appears flipped. To handle
/// that, use the option to flip textures coordinates on either saving or
/// loading. By default texture coordinates are flipped since this seems
/// the convention found on test cases collected on the web.
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// In the high level interface, shapes are indexed meshes and are described
/// by arrays of vertex indices for points/lines/triangles and arrays for vertex
/// positions, normals, texcoords, color and radius. The latter two as
/// extensions. Since OBJ is a complex formats that does not match well with
/// current GPU rendering / path tracing algorithms, we adopt a simplification
/// similar to other single file libraries:
/// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
///   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
///   vertex dusplication happens thought for same triplets
/// 2. we split shapes on changes to groups and materials, instead of keeping
///   per-face group/material data; this makes the data usable right away in
///   a GPU viewer; this is not a major limitation if we accept the previous
///   point that already changes shapes topology.
///
///
/// USAGE FOR WRITING:
///
/// 1. include this file (more compilation options below)
/// 2. [LOW-LEVEL INTERFACE] fill an obj object with your scene data and save
///    the obj/mtl pair with save_obj()
///    ok = save_obj(filename, obj, error message)
/// 3. [HIGH_LEVEL INTERFACE] create a flattened scene object and turn into an
///    obj with unflatten_obj()
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
///
/// Texture loading depends on yocto_image.h.
///
/// If the texture loading dependency is not desired, it can be disabled by
/// defining YOBJ_NO_IMAGE before including this file.
///
///
/// HISTORY:
/// - v 0.12: change texture loading by flipping uvs rather than images
/// - v 0.11: use yocto_image for texture handling.
/// - v 0.10: switch to .h/.cpp pair
/// - v 0.9: bug fixes and optionally texture skipping
/// - v 0.8: high level interface uses grouping
/// - v 0.7: doxygen comments
/// - v 0.6: bug fixes
/// - v 0.5: removed options to force image formats (image library not reliable)
/// - v 0.4: [major API change] move to modern C++ interface
/// - v 0.3: new API internals and C++ interface
/// - v 0.2: removal of C interface
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace yobj {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#ifndef _YOBJ_H_
#define _YOBJ_H_

// compilation options
#ifdef YOBJ_INLINE
#define YOBJ_API inline
#else
#define YOBJ_API
#endif

#include <array>
#include <cmath>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------

namespace yobj {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float16 = std::array<float, 16>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;

/// @name low level API
/// @{

///
/// Face vertex
///
struct vert {
    int pos;       ///< position
    int texcoord;  ///< texcoord
    int norm;      ///< normal
    int color;     ///< color [extension]
    int radius;    ///< radius [extension]

    /// Constructor (copies members initializing missing ones to -1)
    vert(int pos = -1, int texcoord = -1, int norm = -1, int color = -1,
        int radius = -1)
        : pos(pos)
        , texcoord(texcoord)
        , norm(norm)
        , color(color)
        , radius(radius) {}
};

///
/// Elelemnt vertex indices
///
struct elem {
    /// element type
    enum struct etype : uint16_t {
        point = 1,  ///< lists of points
        line = 2,   ///< polylines
        face = 3    ///< polygon faces
    };

    uint32_t start;  ///< starting vertex index
    etype type;      ///< element type
    uint16_t size;   ///< number of vertices

#ifdef _WIN32
    elem(uint32_t start_, etype type_, uint16_t size_)
        : start(start_), type(type_), size(size_) {}
#endif
};

///
/// Element group
///
struct elem_group {
    // group data
    std::string matname;    ///< material name
    std::string groupname;  ///< group name

    // element data
    std::vector<vert> verts;  ///< element vertices
    std::vector<elem> elems;  ///< element faces
};

///
/// Obj object
///
struct object {
    // object data
    std::string name;  ///< object name

    // element data
    std::vector<elem_group> elems;  ///< element groups
};

///
/// OBJ material
///
struct material {
    // whole material data
    std::string name;  ///< material name
    int illum = 0;     ///< MTL illum mode

    // color information
    float3 ke = {0, 0, 0};  ///< emission color
    float3 ka = {0, 0, 0};  ///< ambient color
    float3 kd = {0, 0, 0};  ///< diffuse color
    float3 ks = {0, 0, 0};  ///< specular color
    float3 kr = {0, 0, 0};  ///< reflection color
    float3 kt = {0, 0, 0};  ///< transmision color
    float ns = 1;           ///< phong exponent for ks
    float ior = 1;          ///< index of refraction
    float op = 1;           ///< opacity

    // texture names for the above properties
    std::string ke_txt;    ///< emission texture
    std::string ka_txt;    ///< ambient texture
    std::string kd_txt;    ///< diffuse texture
    std::string ks_txt;    ///< specular texture
    std::string kr_txt;    ///< reflection texture
    std::string kt_txt;    ///< transmission texture
    std::string ns_txt;    ///< specular exponent texture
    std::string op_txt;    ///< opacity texture
    std::string ior_txt;   ///< index of refraction
    std::string bump_txt;  ///< bump map texture (heighfield)
    std::string disp_txt;  ///< displacement map texture (heighfield)
    std::string norm_txt;  ///< normal map texture
};

///
/// Camera [extension]
///
struct camera {
    std::string name;  ///< camera name
    float16 xform = {
        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};  ///< camera transform
    bool ortho = false;                ///< orthografic camera
    float yfov = 2 * std::atan(0.5f);  ///< vertical field of view
    float aspect = 16.0f / 9.0f;       ///< aspect ratio
    float aperture = 0;                ///< lens aperture
    float focus = 1;                   ///< focus distance
};

///
/// Environment [extension]
///
struct environment {
    std::string name;  ///< environment name
    float16 xform = {
        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};  /// transform
    std::string matname;                                  /// material name
};

///
/// OBJ asset
///
struct obj {
    // vertex data
    std::vector<float3> pos;       /// vertex positions
    std::vector<float3> norm;      /// vertex normals
    std::vector<float2> texcoord;  /// vertex texcoord
    std::vector<float3> color;     /// vertex color [extension]
    std::vector<float> radius;     /// vertex radius [extension]

    // scene objects
    std::vector<object> objects;            /// objects
    std::vector<material> materials;        /// materials
    std::vector<camera> cameras;            /// cameras [extension]
    std::vector<environment> environments;  /// env maps [extension]
};

///
/// IO Exception.
///
struct obj_exception : std::exception {
    /// constructor with error message
    obj_exception(const std::string& errmsg) : _errmsg(errmsg) {}

    /// retieval of error message
    virtual const char* what() const throw() { return _errmsg.c_str(); }

   private:
    std::string _errmsg;
};

///
/// Load OBJ
///
/// Parameters:
/// - filename: filename
/// - flip_texcoord: whether to flip the v coordinate
///
/// Return:
/// - obj
///
/// Throws:
/// - io_exception: read/write exception
///
YOBJ_API obj* load_obj(const std::string& filename, bool flip_texcoord = true);

///
/// Load MTL
///
/// Parameters:
/// - filename: filename
///
/// Return:
/// - loaded materials
///
/// Throws:
/// - io_exception: read/write exception
///
YOBJ_API std::vector<material> load_mtl(const std::string& filename);

///
/// Save OBJ
///
/// Parameters:
/// - filename: filename
/// - asset: obj data to save
/// - flip_texcoord: whether to flip the v coordinate
///
/// Throws:
/// - io_exception: read/write exception
///
YOBJ_API void save_obj(
    const std::string& filename, const obj* asset, bool flip_texcoord = true);

///
/// Save MTL (@deprecated interface)
///
/// Throws:
/// - io_exception: read/write exception
///
YOBJ_API void save_mtl(
    const std::string& filename, const std::vector<material>& materials);

/// @}

// -----------------------------------------------------------------------------
// HIGH-LEVEL INTERFACE
// -----------------------------------------------------------------------------

/// @name high level API
/// @{

///
/// Mesh primitives. May contain only one of the points/lines/triangles.
///
struct fl_primitives {
    /// name of the group that enclosed it
    std::string name = "";
    /// material id (-1 if not found)
    int material = -1;

    // shape elements
    std::vector<int> points;      ///< points
    std::vector<int2> lines;      ///< lines
    std::vector<int3> triangles;  ///< triangles

    // vertex data
    std::vector<float3> pos;       ///< per-vertex position (3 float)
    std::vector<float3> norm;      ///< per-vertex normals (3 float)
    std::vector<float2> texcoord;  ///< per-vertex texcoord (2 float)
    std::vector<float3> color;     ///< [extension] per-vertex color (3 float)
    std::vector<float> radius;     ///< [extension] per-vertex radius (1 float)
};

///
/// Scene geometry
///
struct fl_mesh {
    // name
    std::string name;

    /// primitives
    std::vector<int> primitives;
};

///
/// Scene Material
///
struct fl_material {
    // whole material data
    std::string name;  ///< material name

    // color information
    float3 ke = {0, 0, 0};  ///< emission color
    float3 kd = {0, 0, 0};  ///< diffuse color
    float3 ks = {0, 0, 0};  ///< specular color
    float3 kt = {0, 0, 0};  ///< transmission color
    float rs = 0.0001;      ///< roughness

    // indices in the texture array (-1 if not found)
    int ke_txt = -1;    ///< emission texture index
    int kd_txt = -1;    ///< diffuse texture index
    int ks_txt = -1;    ///< specular texture index
    int kt_txt = -1;    ///< transmission texture index
    int rs_txt = -1;    ///< roughness texture index
    int norm_txt = -1;  ///< normal texture index
};

///
/// Scene Texture
///
struct fl_texture {
    std::string path;                  ///< path
    int width = 0, height = 0;         ///< if loaded, image width and hieght
    int ncomp = 0;                     ///< if loaded, number of component (1-4)
    std::vector<float> dataf;          ///< if loaded, pixel data for HDRs
    std::vector<unsigned char> datab;  ///< if loaded, pixel data for LDRs
};

///
/// Scene Camera
///
struct fl_camera {
    std::string name;  ///< name
    float16 xform = {
        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};  ///< transform
    bool ortho = false;                                   ///< ortho cam
    float yfov = 2;               ///< vertical field of view
    float aspect = 16.0f / 9.0f;  ///< aspect ratio
    float aperture = 0;           ///< lens aperture
    float focus = 1;              ///< focus distance
};

///
/// Envinonment map
///
struct fl_environment {
    std::string name;  ///< name
    int matid = -1;  ///< index of material in material array (-1 if not found)
    float16 xform = {
        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};  ///< transform
};

///
/// Scene
///
struct fl_obj {
    std::vector<fl_primitives*> primitives;     ///< shape primitives
    std::vector<fl_mesh*> meshes;               ///< mesh array
    std::vector<fl_material*> materials;        ///< material array
    std::vector<fl_texture*> textures;          ///< texture array
    std::vector<fl_camera*> cameras;            ///< camera array
    std::vector<fl_environment*> environments;  ///< environment array

    ~fl_obj();
};

///
/// Load an asset
///
/// Parameters:
/// - obj: obj to be flattened
///
/// Returns:
/// - flattened scene
///
YOBJ_API fl_obj* flatten_obj(const obj* asset);

///
/// Save an asset
///
/// Parameters:
/// - scene: scene to unflatten
///
/// Returns:
/// - obj
///
YOBJ_API obj* unflatten_obj(const fl_obj* scene);

///
/// Loads textures for an scene.
///
/// Parameters:
/// - scene: scene to load textures into
/// - dirname: base directory name for texture files
/// - skip_missing: whether to skip missing textures or throw an expection
///
/// Throws:
/// - obj_exception
///
YOBJ_API void load_textures(
    fl_obj* scene, const std::string& dirname, bool skip_missing = false);

/// @}

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YOBJ_INLINE
#include "yocto_obj.cpp"
#endif

#endif
