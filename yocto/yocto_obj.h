//
// YOCTO_OBJ: Wavefront OBJ/MTL loader and writer with support for points,
// lines, triangles and general polygons and all materials properties.
// Contains also a few extension to eqasily create demos such as per-vertex
// color and radius, cameras and envmaps.
//

//
// USAGE FOR READING:
//
// 0. include this file (more compilation options below)
// 1. load an obj with yo_load_obj
//   - loads an obj from disk including its associate mtl files
//   - returns a parsed scene data structure described below
//   - optionally support triangulation on loads that makes the loader
//     faster (use it alwayd if you do not need quads/polys/polylines)
//   - extension can be optionally enabled
//   scene = yo_load_obj(filename, triangulate, enable_extensions)
//   1.a. optionally load textures data as float arrays with
//   yo_load_textures(scene, scene_filename, req_comp)
// 2. access the data directly from the returned scene object
//   - has five main arrays: shapes, materials, textures, cameras, envmaps
//   e.g. for(int i = 0; i < scene->nshapes; i ++) scene->shapes[i].XXX
// 3. cleanup with yo_free_scene
//   - you have to do this for each shape bvh and the scene bvh
//   yb_free_bvh(bvh)
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles, quads),
// and arrays for vertex positions, normals, texcoords, color and radius.
// The latter two as extensions.
//
// Quad meshes are experimental and might go away in future realeases. If
// you can, please use triangles. Quads are treated as two triangles (v0,v1,v3)
// and (v2,v3,v1). Quads with v2 == v3 are degenerate and represent one
// triangle, thus quads meshes can also represent mixtures of triangle and
// quads. This follows Intel's Embree API.
//
// Faces in the scene have the same number of elements for points (1),
// lines (2), triangles (3) and quads (4 with the above note). We also
// support general polygons and polines with arbitrary number of faces.
// To avdoi wasting memory, these are saved sequentially where the first int
// in the element is the number od vertices. While this does not allow
// random access, it saves significant memory and avoid pointer chasing.
//
// Since OBJ is a complex formats that does not match well with current
// GPU rendering / path tracing algorithms, we adopt a simplification similar
// to other single file libraries.
//
// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
//   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
//   vertex dusplication happens thought for same triplets
// 2. we split shapes on changes to groups and materials, instead of keeping
//   per-face group/material data; this makes the data usable right away in
//   a GPU viewer; this is not a major limitation if we accept the previous
//   point that already changes shapes topology.
//

//
// USAGE FOR WRITING:
//
// 0. include this file (more compilation options below)
// 1. fill a yo_scene data with your scene data
//    - note that if your shape data is layed out in memore as ours,
//    - then no copy is needed, just set the pointers
// 2. save the obj/mtl pair with yo_load_obj
//    yo_save_obj(filename, scene, enable_extensions)
// 3. If you copies memory over, clear it with yo_free_scene
//

//
// USAGE FOR BINARY DUMPS:
//
// 1. you can also have binary dumps used for fast data access with
//    scene = yo_load_objbin(filename, enable_extensions) and
//    yo_save_objbin(filename, scene, ext)
// 2. These files are just binary dumps, so should not be used for
//    archival but as a speed up to avoid ASCII serializatiion/deserialization
//

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use.
// To use the library as a .h/.c pair do the following:
// - to use as a .h, just #define YO_NOINLINE before including this file
// - to use as a .c, just #define YO_IMPLEMENTATION before including this file
// To disable texture loading code, #define YO_NOIMG (uses stb_image.h to load).
// Internal yocto_obj uses a fixed size hash table to resolve unique vertices
// in OBJ. Use YO_HASHSIZE to change the number of bucks.
//

//
// HISTORY:
// - v 0.1: handles larger files but allocates more up front memory
// - v 0.0: initial release

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YO_H_
#define _YO_H_

#ifndef YO_NOINLINE
#define YO_API static inline
#else
#ifdef __cplusplus
#define YO_API extern "C"
#else
#define YO_API
#endif
#endif

#include <stdbool.h>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Types of geometric primitives
//
enum {
    yo_etype_null = 0,       // invalid prim to indicate parsing erros
    yo_etype_point = 1,      // points
    yo_etype_line = 2,       // lines
    yo_etype_triangle = 3,   // triangles
    yo_etype_quad = 4,       // quads
    yo_etype_polyline = 12,  // polylines
    yo_etype_polygon = 13    // polygons
};

//
// Geometric shape
//
typedef struct yo_shape {
    // whole shape data
    char* name;       // shape name
    char* groupname;  // groupname (unique group for each shape object)
    char* matname;    // material name
    int matid;        // index in the material array (-1 if not found)

    // shape elements
    int nelems;  // number of elements (point, lines, triangles, etc.)
    int* elem;   // per-element vertex indices
    int etype;   // element type from the above enum

    // vertex data
    int nverts;       // number of vertices
    float* pos;       // per-vertex position (3 float)
    float* norm;      // per-vertex normals (3 float)
    float* texcoord;  // per-vertex texcoord (2 float)
    float* color;     // [extension] per-vertex color (3 float)
    float* radius;    // [extension] per-vertex radius (1 float)
    float* xform;     // [extension] 4x4 transform matrix
} yo_shape;

//
// Material
//
typedef struct yo_material {
    // whole material data
    char* name;  // material name
    int illum;   // MTL illum mode

    // color information
    float ke[3];  // emission color
    float ka[3];  // ambient color
    float kd[3];  // diffuse color
    float ks[3];  // specular color
    float kr[3];  // reflection color
    float kt[3];  // transmision color
    float ns;     // phong exponent for ks
    float ior;    // index of refraction
    float op;     // opacity

    // texture names for the above properties
    char* ke_txt;
    char* ka_txt;
    char* kd_txt;
    char* ks_txt;
    char* kr_txt;
    char* kt_txt;
    char* ns_txt;
    char* op_txt;
    char* ior_txt;
    char* bump_txt;  // bump map texture (heighfield)
    char* disp_txt;  // displacement map texture (heighfield)

    // indices in the texture array (-1 if not found)
    int ke_txtid;
    int ka_txtid;
    int kd_txtid;
    int ks_txtid;
    int kr_txtid;
    int kt_txtid;
    int ns_txtid;
    int op_txtid;
    int ior_txtid;
    int bump_txtid;
    int disp_txtid;
} yo_material;

//
// [Extension] Texture
//
typedef struct yo_texture {
    char* path;         // path
    int width, height;  // if loaded, image width and hieght
    int ncomp;          // if loaded, number of component (1-4)
    float* pixels;      // if loaded, pixel data
} yo_texture;

//
// [Extension] Camera represented as a lookat.
//
typedef struct yo_camera {
    char* name;           // name
    float from[3];        // camera position
    float to[3];          // camera focus location
    float up[3];          // camera up vector
    float width, height;  // image plane width and height
    float aperture;       // lens aperture
} yo_camera;

//
// [Extension] Envinonment map in latlong format
//
typedef struct yo_env {
    char* name;     // name
    char* matname;  // material name (where only ke, ke_txt are valid)
    int matid;      // index of material in material array (-1 if not found)
    float from[3], to[3], up[3];  // lookat transform data as in yo_camera
} yo_env;

//
// Scene
//
typedef struct yo_scene {
    int nshapes;             // number of shapes
    yo_shape* shapes;        // shape array
    int nmaterials;          // number of materials
    yo_material* materials;  // material array
    int ntextures;           // number of textures
    yo_texture* textures;    // texture array
    int ncameras;            // number of cameras
    yo_camera* cameras;      // camera array
    int nenvs;               // number of environment
    yo_env* envs;            // environment array
} yo_scene;

//
// Loads a scene from disk
//
// Parameters:
// - filename: scene filename
// - truangulate: whether to triagulate on load (fan-style)
// - ext: enable extensions
//
// Returns:
// - loaded scene or NULL for error
//
YO_API yo_scene*
yo_load_obj(const char* filename, bool triangulate, bool ext);

//
// Loads a binary scene dump from disk
//
// Parameters:
// - filename: scene filename
// - ext: enable extensions
//
// Returns:
// - loaded scene or NULL for error
//
YO_API yo_scene*
yo_load_objbin(const char* filename, bool ext);

//
// Saves a scene to disk
//
// Parameters:
// - filename: scene filename
// - scene: scene to save
// - ext: enable extensions
//
// Returns:
// - true if ok
//
YO_API bool
yo_save_obj(const char* filename, const yo_scene* scene, bool ext);

//
// Saves a binary scene dump to disk
//
// Parameters:
// - filename: scene filename
// - scene: scene to save
// - ext: enable extensions
//
// Returns:
// - true if ok
//
YO_API bool
yo_save_objbin(const char* filename, const yo_scene* scene, bool ext);

//
// Free scene data.
//
YO_API void
yo_free_scene(yo_scene* scene);

//
// Loads textures.
//
// Parameters:
// - scene: scene to load into
// - filename: scene filename, used to resolve path references
// - req_comp: 0 for default or 1-4 to force all textures to have the given
//    number of components
//
#ifndef YO_NOIMG
YO_API void
yo_load_textures(yo_scene* scene, const char* filename, int req_comp);
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YO_NOINLINE) || defined(YO_IMPLEMENTATION)

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// LOW-LEVEL SUPPORT FOR FIXED VECTORS AND GROWABLE ARRAYS
// -----------------------------------------------------------------------------

//
// 2d float vectors
//
typedef struct { float x, y; } yo__vec2f;

//
// 3d float vectors
//
typedef struct { float x, y, z; } yo__vec3f;

//
// 4x4 float transform matrices
//
typedef struct { float m[16]; } yo__mat4f;

//
// OBJ vertex reference triplet (pos,texcoord,norm) with extension for
// color and radius indices. Contains also the vertex unique index in the
// flattened array.
//
typedef struct { int pos, texcoord, norm, color, radius, vid; } yo__vert;

//
// Vertex hash table to avoid duplicating vertices.
//
#define yo__vhash_size 1048576
typedef struct yo__vhash {
    int nverts;                   // numner of vertices
    int s[yo__vhash_size];        // size of the hash buckets
    yo__vert* v[yo__vhash_size];  // bucket data (one vertex for
                                  // each unique face index)
} yo__vhash;

//
// Round up to next power of two to use with growable arrays without
// keeping the capicity explicitly.
//
static inline int
yo__vector_capacity(int s) {
    if (!s) return 0;
    // round to next pow2
    // http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    uint32_t n = s;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

// macro to create new objects
#define yo__new(t) ((t*)calloc(1, sizeof(t)))

//
// Grow an array allocated with a capacity that is the next power of two of the
// size.
//
#define yo__vector_grow(t, d, s, n)                                            \
    {                                                                          \
        if (!s || yo__vector_capacity(s + n) > yo__vector_capacity(s)) {       \
            d = (t*)realloc(d, sizeof(t) * yo__vector_capacity(s + n));        \
            memset(d + s, 0, sizeof(*d) * (yo__vector_capacity(s + n) - s));   \
        }                                                                      \
        s += n;                                                                \
    }

//
// Pushback a value into a groable array
//
#define yo__pushback(t, d, s, vv)                                              \
    {                                                                          \
        if (!(s) || yo__vector_capacity((s) + 1) > yo__vector_capacity((s))) { \
            (d) = (t*)realloc(d, sizeof(t) * yo__vector_capacity((s) + 1));    \
        }                                                                      \
        (d)[(s)] = (vv);                                                       \
        (s)++;                                                                 \
    }

//
// Trim growable array to the right size
//
#define yo__trim(t, v, n)                                                      \
    {                                                                          \
        if (v) v = (t*)realloc((v), sizeof(t) * (n));                          \
    }

//
// strdup dropin replacement that handle gracefully NULL strings
//
static inline char*
yo__strdup(const char* str) {
    if (!str) return 0;
    char* ret = (char*)calloc(strlen(str) + 1, sizeof(char));
    strcpy(ret, str);
    return ret;
}

//
// String comparison without case.
// modified from http://clc-wiki.net/wiki/C_standard_library%3astring.h%3astrcmp
//
static inline int
yo__stricmp(const char* s1, const char* s2) {
    while (*s1 && (tolower(*s1) == tolower(*s2))) s1++, s2++;
    return tolower(*(const unsigned char*)s1) -
           tolower(*(const unsigned char*)s2);
}

// -----------------------------------------------------------------------------
// OBJ LOADING
// -----------------------------------------------------------------------------

//
// Free scene memory.
//
YO_API void
yo_free_scene(yo_scene* scene) {
    for (int i = 0; i < scene->ncameras; i++) {
        yo_camera* cam = &scene->cameras[i];
        if (cam->name) free(cam->name);
    }
    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = &scene->shapes[i];
        if (shape->name) free(shape->name);
        if (shape->matname) free(shape->matname);
        if (shape->elem) free(shape->elem);
        if (shape->pos) free(shape->pos);
        if (shape->norm) free(shape->norm);
        if (shape->texcoord) free(shape->texcoord);
    }
    for (int i = 0; i < scene->nmaterials; i++) {
        yo_material* mat = &scene->materials[i];
        if (mat->name) free(mat->name);
        char* txts[10] = { mat->ke_txt,  mat->ka_txt, mat->kd_txt,
                           mat->ks_txt,  mat->kr_txt, mat->kt_txt,
                           mat->ns_txt,  mat->op_txt, mat->bump_txt,
                           mat->disp_txt };
        for (int t = 0; t < 10; t++)
            if (txts[t]) free(txts[t]);
    }
    for (int i = 0; i < scene->nenvs; i++) {
        yo_env* env = &scene->envs[i];
        if (env->name) free(env->name);
        if (env->matname) free(env->matname);
    }
    for (int i = 0; i < scene->ntextures; i++) {
        yo_texture* txt = &scene->textures[i];
        if (txt->path) free(txt->path);
        if (txt->pixels) free(txt->pixels);
    }
}

//
// During parsing, flashes a shape into the scene if elements are present.
//
YO_API yo_shape*
yo__flush_shape(yo_scene* scene, yo__vhash* vhash) {
    // exit if nothing to do
    yo_shape* shape = scene->shapes + scene->nshapes - 1;
    if (!shape->nelems) return shape;

    // trim vertices
    yo__trim(float, shape->pos, shape->nverts * 3);
    yo__trim(float, shape->norm, shape->nverts * 3);
    yo__trim(float, shape->texcoord, shape->nverts * 2);
    yo__trim(float, shape->color, shape->nverts * 3);
    yo__trim(float, shape->radius, shape->nverts);

    // handle simple cases for elements
    if (shape->etype == yo_etype_point || shape->etype == yo_etype_line ||
        shape->etype == yo_etype_triangle || shape->etype == yo_etype_quad) {
        shape->nelems /= shape->etype;
        yo__trim(int, shape->elem, shape->nelems * shape->etype);
    } else if (shape->etype == yo_etype_polygon ||
               shape->etype == yo_etype_polyline) {
        // tries to compress generic polygon and polylines
        // find size
        int nelems = shape->nelems;
        int* elem = shape->elem;
        shape->nelems = 0;
        int maxf = INT_MIN, minf = INT_MAX;
        for (int f = 0; f < nelems;) {
            int nf = elem[f];
            if (nf > maxf) maxf = nf;
            if (nf < minf) minf = nf;
            f += nf + 1;
            shape->nelems++;
        }
        assert(minf > 0);

        // make triangle and quad meshes
        if (minf == maxf && maxf < 5) {
            shape->etype = maxf;
            shape->elem =
                (int*)calloc(shape->nelems * shape->etype, sizeof(int));
            for (int e = 0; e < shape->nelems; e++)
                memcpy(shape->elem + e * maxf, elem + e * (maxf + 1) + 1,
                       sizeof(int) * maxf);
            free(elem);
        } else if (minf == 3 && maxf == 4) {
            shape->etype = maxf;
            shape->elem =
                (int*)calloc(shape->nelems * shape->etype, sizeof(int));
            int epos = 0;
            for (int f = 0; f < nelems;) {
                int nf = elem[f];
                memcpy(shape->elem + epos, elem + f + 1, sizeof(int) * nf);
                if (nf == 3)
                    *(shape->elem + epos + 3) = *(shape->elem + epos + 2);
                f += nf + 1;
                epos += 4;
            }
            free(elem);
        } else {
            yo__trim(int, shape->elem, nelems);
        }
    } else {
        assert(false);
    }

    // clear buffers
    vhash->nverts = 0;
    for (int i = 0; i < yo__vhash_size; i++) vhash->s[i] = 0;

    // get a new shape
    yo__vector_grow(yo_shape, scene->shapes, scene->nshapes, 1);
    return scene->shapes + scene->nshapes - 1;
}

//
// Splits a string into an array of strings on whitespace with Python split
// semantic. Modifies original string to avoid allocation.
//
static inline int
yo__splitws(char* str, char** splits, int maxsplits) {
    int n = 0;
    while (*str && n < maxsplits) {
        if (isspace(*str)) {
            *str = 0;
        } else {
            if (n == 0 || !(*(str - 1))) {
                splits[n] = str;
                n++;
            }
        }
        str++;
    }
    return n;
}

//
// Add an empty material.
//
static inline void
yo__add_empty_material(yo_scene* scene, const char* name) {
    yo__vector_grow(yo_material, scene->materials, scene->nmaterials, 1);
    yo_material* mat = scene->materials + scene->nmaterials - 1;
    mat->name = yo__strdup(name);
}

//
// Add a camera from OBJ vertices.
//
static inline void
yo__add_camera(yo_scene* scene, const char* name, yo__vert from, yo__vert to,
               yo__vec3f* pos, yo__vec3f* norm, yo__vec2f* texcoord,
               yo__vhash* vhash) {
    yo__vector_grow(yo_camera, scene->cameras, scene->ncameras, 1);
    yo_camera* cam = scene->cameras + scene->ncameras - 1;

    cam->name = yo__strdup(name);
    *(yo__vec3f*)cam->from = pos[from.pos];
    *(yo__vec3f*)cam->to = pos[to.pos];
    *(yo__vec3f*)cam->up =
        (from.norm >= 0) ? norm[from.norm] : (yo__vec3f){ 0, 1, 0 };
    cam->width = (to.texcoord >= 0) ? texcoord[to.texcoord].x : 1;
    cam->height = (to.texcoord >= 0) ? texcoord[to.texcoord].y : 1;
    cam->aperture = (from.texcoord >= 0) ? texcoord[from.texcoord].x : 0;

    // clear buffers
    vhash->nverts = 0;
    for (int i = 0; i < yo__vhash_size; i++) vhash->s[i] = 0;
}

//
// Add an environment map from OBJ vertices.
//
static inline void
yo__add_env(yo_scene* scene, const char* name, const char* matname,
            yo__vert from, yo__vert to, yo__vec3f* pos, yo__vec3f* norm,
            yo__vhash* vhash) {
    yo__vector_grow(yo_env, scene->envs, scene->nenvs, 1);
    yo_env* env = scene->envs + scene->nenvs - 1;

    env->name = yo__strdup(name);
    env->matname = yo__strdup(matname);
    *(yo__vec3f*)env->from = pos[from.pos];
    *(yo__vec3f*)env->to = pos[to.pos];
    *(yo__vec3f*)env->up =
        (from.norm >= 0) ? norm[from.norm] : (yo__vec3f){ 0, 1, 0 };

    // clear buffers
    vhash->nverts = 0;
    for (int i = 0; i < yo__vhash_size; i++) vhash->s[i] = 0;
}

// parses one float
static inline float
yo__parse_float(char** tok) {
    return atof(tok[0]);
}

// parses two floats
static inline yo__vec2f
yo__parse_float2(char** tok) {
    return (yo__vec2f){ (float)atof(tok[0]), (float)atof(tok[1]) };
}

// parses three floats
static inline yo__vec3f
yo__parse_float3(char** tok) {
    return (yo__vec3f){ (float)atof(tok[0]), (float)atof(tok[1]),
                        (float)atof(tok[2]) };
}

// parses 16 floats
static inline yo__mat4f
yo__parse_mat4f(char** tok) {
    yo__mat4f m = { { 0 } };
    for (int i = 0; i < 16; i++) m.m[i] = (float)atof(tok[i]);
    return m;
}

// parses an OBJ vertex triplet (or quintuplet with extensions); handle
// nagative indices directly
static inline yo__vert
yo__parse_vert(char* str, yo__vhash* vhash, yo__vert vl) {
    // parse triplet
    char* splits[] = { str, 0, 0, 0, 0 };
    int ns = 1;
    yo__vert v = { -1, -1, -1, -1, -1 };
    while (*str) {
        if (*str == '/') {
            *str = 0;
            if (ns < 5) splits[ns++] = str + 1;
        }
        str++;
    }
    int* f = &v.pos;
    int* l = &vl.pos;
    for (int i = 0; i < 5; i++) {
        if (!splits[i]) {
            f[i] = -1;
            continue;
        }
        f[i] = (int)atoi(splits[i]);
        f[i] = (f[i] < 0) ? l[i] + f[i] : f[i] - 1;
    }

    // determine position vid using vertex hash
    int pos = -1;
    int hidx = v.pos % yo__vhash_size;
    for (int i = 0; i < vhash->s[hidx] && pos < 0; i++) {
        if (v.pos == vhash->v[hidx][i].pos &&
            v.texcoord == vhash->v[hidx][i].texcoord &&
            v.norm == vhash->v[hidx][i].norm &&
            v.color == vhash->v[hidx][i].color &&
            v.radius == vhash->v[hidx][i].radius)
            pos = i;
    }

    // found, can exit
    if (pos >= 0) return vhash->v[hidx][pos];

    // insert in vhash
    v.vid = vhash->nverts;
    yo__pushback(yo__vert, vhash->v[hidx], vhash->s[hidx], v);
    vhash->nverts++;

    return v;
}

// add a unique vertex to a parsed shape
static inline void
yo__add_shape_vert(yo_shape* shape, yo__vert v, yo__vec3f* pos, yo__vec3f* norm,
                   yo__vec2f* texcoord, yo__vec3f* color, float* radius) {
    // check already added
    if (v.vid < shape->nverts) return;
    // TODO: assert for malformed stuff
    if (v.pos >= 0) {
        int npos = shape->nverts * 3;
        yo__vector_grow(float, shape->pos, npos, 3);
        shape->pos[shape->nverts * 3 + 0] = pos[v.pos].x;
        shape->pos[shape->nverts * 3 + 1] = pos[v.pos].y;
        shape->pos[shape->nverts * 3 + 2] = pos[v.pos].z;
    }
    if (v.norm >= 0) {
        int nnorm = shape->nverts * 3;
        yo__vector_grow(float, shape->norm, nnorm, 3);
        shape->norm[shape->nverts * 3 + 0] = norm[v.norm].x;
        shape->norm[shape->nverts * 3 + 1] = norm[v.norm].y;
        shape->norm[shape->nverts * 3 + 2] = norm[v.norm].z;
    }
    if (v.texcoord >= 0) {
        int ntexcoord = shape->nverts * 2;
        yo__vector_grow(float, shape->texcoord, ntexcoord, 2);
        shape->texcoord[shape->nverts * 2 + 0] = texcoord[v.texcoord].x;
        shape->texcoord[shape->nverts * 2 + 1] = texcoord[v.texcoord].y;
    }
    if (v.color >= 0) {
        int ncolor = shape->nverts * 3;
        yo__vector_grow(float, shape->color, ncolor, 3);
        shape->color[shape->nverts * 3 + 0] = color[v.color].x;
        shape->color[shape->nverts * 3 + 1] = color[v.color].y;
        shape->color[shape->nverts * 3 + 2] = color[v.color].z;
    }
    if (v.radius >= 0) {
        int nradius = shape->nverts;
        yo__vector_grow(float, shape->radius, nradius, 1);
        shape->radius[shape->nverts] = radius[v.radius];
    }
    shape->nverts += 1;
}

//
// loads an MTL file
//
static inline yo_scene*
yo__load_mtl(yo_scene* scene, const char* filename) {
    FILE* mfile = fopen(filename, "rt");
    if (!mfile) return 0;

    char mline[4096];
    char* mtok[10];
    int mlinenum = 0;

    yo_material* mat = 0;

    // foreach line, splits the line by whitespaces and parses the data
    // directly in the material
    while (fgets(mline, 4096, mfile)) {
        mlinenum += 1;
        int mntok = yo__splitws(mline, mtok, 10);

        if (!mntok) {
            continue;
        } else if (mtok[0][0] == '#' || mtok[0][0] == '/') {
            continue;
        } else if (!strcmp(mtok[0], "newmtl")) {
            yo__add_empty_material(scene, mtok[1]);
            mat = &scene->materials[scene->nmaterials - 1];
        } else if (!strcmp(mtok[0], "illum")) {
            mat->illum = atoi(mtok[1]);
        } else if (!strcmp(mtok[0], "Ke")) {
            *(yo__vec3f*)mat->ke = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Ka")) {
            *(yo__vec3f*)mat->ka = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Kd")) {
            *(yo__vec3f*)mat->kd = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Ks")) {
            *(yo__vec3f*)mat->ks = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Kr")) {
            *(yo__vec3f*)mat->kr = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Tr")) {
            *(yo__vec3f*)mat->kt = yo__parse_float3(mtok + 1);
        } else if (!strcmp(mtok[0], "Ns")) {
            mat->ns = yo__parse_float(mtok + 1);
        } else if (!strcmp(mtok[0], "d")) {
            mat->op = yo__parse_float(mtok + 1);
        } else if (!strcmp(mtok[0], "Tr")) {
            mat->op = yo__parse_float(mtok + 1);
        } else if (!strcmp(mtok[0], "Ni")) {
            mat->ior = yo__parse_float(mtok + 1);
        } else if (!strcmp(mtok[0], "map_Ke")) {
            mat->ke_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Ka")) {
            mat->ka_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Kd")) {
            mat->kd_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Ks")) {
            mat->ks_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Kr")) {
            mat->kr_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Tr")) {
            mat->kt_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Ns")) {
            mat->ns_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_d")) {
            mat->op_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Tr")) {
            mat->op_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_Ni")) {
            mat->ior_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_bump")) {
            mat->bump_txt = yo__strdup(mtok[1]);
        } else if (!strcmp(mtok[0], "map_disp")) {
            mat->disp_txt = yo__strdup(mtok[1]);
        } else {
            // printf("ignoring value for %s\n", mtok[0]);
        }
    }

    fclose(mfile);

    return scene;
}

//
// Splits a path into component to get directory name
//
static inline void
yo__split_path(const char* filename, char* dirname, char* basename, char* ext) {
    // walk till end keeping the position of '/', '\\' and '.'
    const char *path_sep = 0, *ext_sep = 0;
    for (const char* p = filename; *p; p++) {
        if (*p == '/' || *p == '\\') path_sep = p;
        if (*p == '.') ext_sep = p;
    }

    // copy strings
    if (dirname) {
        if (path_sep) {
            strncpy(dirname, filename, 1 + path_sep - filename);
            if (path_sep) dirname[1 + path_sep - filename] = 0;
        } else
            strcpy(dirname, "");
    }
    if (basename) {
        const char* start = (path_sep) ? path_sep + 1 : filename;
        if (ext_sep) {
            strncpy(basename, start, ext_sep - filename);
            if (ext_sep) basename[ext_sep - start] = 0;
        } else
            strcpy(basename, start);
    }
    if (ext) {
        if (ext_sep)
            strcpy(ext, ext_sep);
        else
            strcpy(ext, "");
    }
}

//
// Sets material indices at the end of file loading
//
static inline void
yo__set_matids(yo_scene* scene) {
    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = scene->shapes + i;
        shape->matid = -1;
        if (shape->matname) {
            for (int j = 0; j < scene->nmaterials && shape->matid < 0; j++) {
                if (!yo__stricmp(scene->materials[j].name, shape->matname))
                    shape->matid = j;
            }
        }
    }
    for (int i = 0; i < scene->nenvs; i++) {
        yo_env* env = scene->envs + i;
        env->matid = -1;
        if (env->matname) {
            for (int j = 0; j < scene->nmaterials && env->matid < 0; j++) {
                if (!yo__stricmp(scene->materials[j].name, env->matname))
                    env->matid = j;
            }
        }
    }
}

//
// Set texture indices at the end of file loading
//
static inline void
yo__set_textures(yo_scene* scene) {
    scene->ntextures = 0;
    scene->textures = 0;
    for (int i = 0; i < scene->nmaterials; i++) {
        yo_material* mat = scene->materials + i;
        char* paths[12] = { mat->ke_txt,  mat->ka_txt,   mat->kd_txt,
                            mat->ks_txt,  mat->kr_txt,   mat->kt_txt,
                            mat->ns_txt,  mat->op_txt,   mat->op_txt,
                            mat->ior_txt, mat->bump_txt, mat->disp_txt };
        int* ids[12] = { &mat->ke_txtid,  &mat->ka_txtid,   &mat->kd_txtid,
                         &mat->ks_txtid,  &mat->kr_txtid,   &mat->kt_txtid,
                         &mat->ns_txtid,  &mat->op_txtid,   &mat->op_txtid,
                         &mat->ior_txtid, &mat->bump_txtid, &mat->disp_txtid };
        for (int j = 0; j < 12; j++) {
            if (!paths[j]) {
                *(ids[j]) = -1;
            } else {
                int pos = -1;
                for (int k = 0; k < scene->ntextures && pos < 0; k++) {
                    if (!strcmp(paths[j], scene->textures[k].path)) pos = k;
                }

                if (pos < 0) {
                    yo__vector_grow(yo_texture, scene->textures,
                                    scene->ntextures, 1);
                    yo_texture* txt = scene->textures + scene->ntextures - 1;
                    txt->path = yo__strdup(paths[j]);
                    pos = scene->ntextures - 1;
                }

                *(ids[j]) = pos;
            }
        }
    }
}

//
// Loads an OBJ file
//
YO_API yo_scene*
yo_load_obj(const char* filename, bool triangulate, bool ext) {
    // prepare scene
    yo_scene* scene = yo__new(yo_scene);

    // vertex scene
    yo__vert obj_nvert = { 0, 0, 0 };
    yo__vec3f* obj_pos = 0;
    yo__vec3f* obj_norm = 0;
    yo__vec2f* obj_texcoord = 0;
    yo__vec3f* obj_color = 0;
    float* obj_radius = 0;

    // current shape scene
    yo__vhash* vhash = (yo__vhash*)calloc(1, sizeof(yo__vhash));
    yo__vector_grow(yo_shape, scene->shapes, scene->nshapes, 1);
    yo_shape* shape = scene->shapes + scene->nshapes - 1;

    // start
    FILE* file = fopen(filename, "rt");
    if (!file) return 0;

    // foreach line, splits the line by whitespaces and parses the data
    // directly in the current shape, emitting shapes when either name,
    // material name, group name or shape element type changes
    char line[4096];
    char* tok[1024];
    int linenum = 0;
    while (fgets(line, 4096, file)) {
        linenum += 1;
        int ntok = yo__splitws(line, tok, 1024);

        if (!ntok) {
            continue;
        } else if (tok[0][0] == '#') {
            continue;
        } else if (!strcmp(tok[0], "v")) {
            yo__vec3f v = yo__parse_float3(tok + 1);
            yo__pushback(yo__vec3f, obj_pos, obj_nvert.pos, v);
        } else if (!strcmp(tok[0], "vt")) {
            yo__vec2f v = yo__parse_float2(tok + 1);
            yo__pushback(yo__vec2f, obj_texcoord, obj_nvert.texcoord, v);
        } else if (!strcmp(tok[0], "vn")) {
            yo__vec3f v = yo__parse_float3(tok + 1);
            yo__pushback(yo__vec3f, obj_norm, obj_nvert.norm, v);
        } else if (!strcmp(tok[0], "vc")) {
            if (ext) {
                yo__vec3f v = yo__parse_float3(tok + 1);
                yo__pushback(yo__vec3f, obj_color, obj_nvert.color, v);
            }
        } else if (!strcmp(tok[0], "vr")) {
            if (ext) {
                float v = yo__parse_float(tok + 1);
                yo__pushback(float, obj_radius, obj_nvert.radius, v);
            }
        } else if (!strcmp(tok[0], "xf")) {
            if (ext) {
                yo__mat4f xf = yo__parse_mat4f(tok + 1);
                shape->xform = (float*)calloc(16, sizeof(float));
                *(yo__mat4f*)shape->xform = xf;
            }
        } else if (!strcmp(tok[0], "c")) {
            if (ext) {
                shape = yo__flush_shape(scene, vhash);
                yo__vert from = yo__parse_vert(tok[1], vhash, obj_nvert);
                yo__vert to = yo__parse_vert(tok[2], vhash, obj_nvert);
                yo__add_camera(scene, shape->name, from, to, obj_pos, obj_norm,
                               obj_texcoord, vhash);
                shape->name = 0;
            }
        } else if (!strcmp(tok[0], "e")) {
            if (ext) {
                shape = yo__flush_shape(scene, vhash);
                yo__vert from = yo__parse_vert(tok[1], vhash, obj_nvert);
                yo__vert to = yo__parse_vert(tok[2], vhash, obj_nvert);
                yo__add_env(scene, shape->name, shape->matname, from, to,
                            obj_pos, obj_norm, vhash);
                shape->name = 0;
                shape->matname = 0;
            }
        } else if (!strcmp(tok[0], "f") && !triangulate) {
            if (shape->etype != yo_etype_polygon) {
                shape = yo__flush_shape(scene, vhash);
            }
            shape->etype = yo_etype_polygon;
            yo__pushback(int, shape->elem, shape->nelems, ntok - 1);
            for (int t = 1; t < ntok; t++) {
                yo__vert v = yo__parse_vert(tok[t], vhash, obj_nvert);
                yo__add_shape_vert(shape, v, obj_pos, obj_norm, obj_texcoord,
                                   obj_color, obj_radius);
                yo__pushback(int, shape->elem, shape->nelems, v.vid);
            }
        } else if (!strcmp(tok[0], "f") && triangulate) {
            if (shape->etype != yo_etype_triangle) {
                shape = yo__flush_shape(scene, vhash);
            }
            shape->etype = yo_etype_triangle;
            int vi0 = 0;
            for (int t = 1; t < ntok; t++) {
                yo__vert v = yo__parse_vert(tok[t], vhash, obj_nvert);
                yo__add_shape_vert(shape, v, obj_pos, obj_norm, obj_texcoord,
                                   obj_color, obj_radius);
                if (t == 1) vi0 = v.vid;
                if (t > 3) {
                    int vil = shape->elem[shape->nelems - 1];
                    yo__pushback(int, shape->elem, shape->nelems, vi0);
                    yo__pushback(int, shape->elem, shape->nelems, vil);
                }
                yo__pushback(int, shape->elem, shape->nelems, v.vid);
            }
        } else if (!strcmp(tok[0], "l") && !triangulate) {
            if (shape->etype != yo_etype_polyline) {
                shape = yo__flush_shape(scene, vhash);
            }
            shape->etype = yo_etype_polyline;
            yo__pushback(int, shape->elem, shape->nelems, ntok - 1);
            for (int t = 1; t < ntok; t++) {
                yo__vert v = yo__parse_vert(tok[t], vhash, obj_nvert);
                yo__add_shape_vert(shape, v, obj_pos, obj_norm, obj_texcoord,
                                   obj_color, obj_radius);
                yo__pushback(int, shape->elem, shape->nelems, v.vid);
            }
        } else if (!strcmp(tok[0], "l") && triangulate) {
            if (shape->etype != yo_etype_line) {
                shape = yo__flush_shape(scene, vhash);
            }
            shape->etype = yo_etype_line;
            for (int t = 1; t < ntok; t++) {
                yo__vert v = yo__parse_vert(tok[t], vhash, obj_nvert);
                yo__add_shape_vert(shape, v, obj_pos, obj_norm, obj_texcoord,
                                   obj_color, obj_radius);
                if (t > 2) {
                    int vil = shape->elem[shape->nelems - 1];
                    yo__pushback(int, shape->elem, shape->nelems, vil);
                }
                yo__pushback(int, shape->elem, shape->nelems, v.vid);
            }
        } else if (!strcmp(tok[0], "p")) {
            if (shape->etype != yo_etype_point) {
                shape = yo__flush_shape(scene, vhash);
            }
            shape->etype = yo_etype_point;
            for (int t = 1; t < ntok; t++) {
                yo__vert v = yo__parse_vert(tok[t], vhash, obj_nvert);
                yo__add_shape_vert(shape, v, obj_pos, obj_norm, obj_texcoord,
                                   obj_color, obj_radius);
                yo__pushback(int, shape->elem, shape->nelems, v.vid);
            }
        } else if (!strcmp(tok[0], "o")) {
            shape = yo__flush_shape(scene, vhash);
            shape->name = yo__strdup((ntok > 1) ? tok[1] : 0);
        } else if (!strcmp(tok[0], "g")) {
            shape = yo__flush_shape(scene, vhash);
            shape->groupname = yo__strdup((ntok > 1) ? tok[1] : 0);
        } else if (!strcmp(tok[0], "usemtl")) {
            shape = yo__flush_shape(scene, vhash);
            shape->matname = yo__strdup((ntok > 1) ? tok[1] : 0);
        } else if (!strcmp(tok[0], "mtllib")) {
            char mfilename[4096];
            yo__split_path(filename, mfilename, 0, 0);
            strcat(mfilename, tok[1]);
            if (!yo__load_mtl(scene, mfilename)) return 0;
        } else {
            // TODO: explicit skips
        }
    }

    // flush and cleanup empty shape if necessary
    yo__flush_shape(scene, vhash);
    scene->nshapes -= 1;

    // close file
    fclose(file);

    // fix ids
    yo__set_matids(scene);
    yo__set_textures(scene);

    // clear
    for (int i = 0; i < yo__vhash_size; i++) {
        if (vhash->v[i]) free(vhash->v[i]);
    }
    free(vhash);
    if (obj_pos) free(obj_pos);
    if (obj_norm) free(obj_norm);
    if (obj_texcoord) free(obj_texcoord);
    if (obj_color) free(obj_color);
    if (obj_radius) free(obj_radius);

    // trim data
    yo__trim(yo_camera, scene->cameras, scene->ncameras);
    yo__trim(yo_shape, scene->shapes, scene->nshapes);
    yo__trim(yo_material, scene->materials, scene->nmaterials);
    yo__trim(yo_texture, scene->textures, scene->ntextures);

    return scene;
}

// -----------------------------------------------------------------------------
// OBJ SAVING
// -----------------------------------------------------------------------------

// write one float prepended by a string
static inline void
yo__fwrite_float(FILE* file, const char* str, float v) {
    fprintf(file, "%s %.6g\n", str, v);
}

// write two floats prepended by a string
static inline void
yo__fwrite_float2(FILE* file, const char* str, float v[2]) {
    fprintf(file, "%s %.6g %.6g\n", str, v[0], v[1]);
}

// write three floats prepended by a string
static inline void
yo__fwrite_float3(FILE* file, const char* str, float v[3]) {
    fprintf(file, "%s %.6g %.6g %.6g\n", str, v[0], v[1], v[2]);
}

// write 16 floats prepended by a string
static inline void
yo__fwrite_float16(FILE* file, const char* str, float v[16]) {
    fprintf(file, "%s", str);
    for (int i = 0; i < 16; i++) fprintf(file, " %.6g", v[i]);
    fprintf(file, "\n");
}

// write a string prepended by another if the string is not NULL
static inline void
yo__fwrite_str(FILE* file, const char* str, const char* s) {
    if (s) fprintf(file, "%s %s\n", str, s);
}

//
// save MTL file
//
static inline bool
yo__save_mtl(const char* filename, const yo_scene* scene) {
    // TODO: failure
    FILE* mfile = fopen(filename, "wt");
    if (!mfile) return false;

    // for each material, dump all the values
    for (int mid = 0; mid < scene->nmaterials; mid++) {
        yo_material* mat = &scene->materials[mid];
        fprintf(mfile, "newmtl %s\n", mat->name);
        fprintf(mfile, "  illum %d\n", mat->illum);
        yo__fwrite_float3(mfile, "  Ke", mat->ke);
        yo__fwrite_float3(mfile, "  Kd", mat->kd);
        yo__fwrite_float3(mfile, "  Ks", mat->ks);
        yo__fwrite_float3(mfile, "  Kr", mat->kr);
        yo__fwrite_float3(mfile, "  Kt", mat->kt);
        yo__fwrite_float(mfile, "  Ns", mat->ns);
        yo__fwrite_float(mfile, "  d", mat->op);
        yo__fwrite_float(mfile, "  Ni", mat->ior);
        yo__fwrite_str(mfile, "  map_Ke", mat->ke_txt);
        yo__fwrite_str(mfile, "  map_Kd", mat->kd_txt);
        yo__fwrite_str(mfile, "  map_Ks", mat->ks_txt);
        yo__fwrite_str(mfile, "  map_Kr", mat->kr_txt);
        yo__fwrite_str(mfile, "  map_Kt", mat->kt_txt);
        yo__fwrite_str(mfile, "  map_Ns", mat->ns_txt);
        yo__fwrite_str(mfile, "  map_d", mat->op_txt);
        yo__fwrite_str(mfile, "  map_Ni", mat->ior_txt);
        yo__fwrite_str(mfile, "  map_bump", mat->bump_txt);
        yo__fwrite_str(mfile, "  map_disp", mat->disp_txt);
        fprintf(mfile, "\n");
    }

    fclose(mfile);

    return true;
}

// write an OBJ vertex triplet using only the indices that are active
static inline void
yo__fwrite_objverts(FILE* file, const char* str, int nv, int* vid,
                    yo__vert voffset, int nto_write, yo__vert to_write) {
    fprintf(file, "%s", str);
    for (int v = 0; v < nv; v++) {
        for (int i = 0; i < nto_write; i++) {
            if ((&to_write.pos)[i])
                fprintf(file, "%c%d", ((i == 0) ? ' ' : '/'),
                        (&voffset.pos)[i] + vid[v]);
            else
                fprintf(file, "%c", '/');
        }
    }
    fprintf(file, "\n");
}

//
// save OBJ
//
YO_API bool
yo_save_obj(const char* filename, const yo_scene* scene, bool ext) {
    char dirname[4096], mfilename[4096];
    yo__split_path(filename, dirname, mfilename, 0);
    strcat(mfilename, ".mtl");

    // write material file
    if (scene->nmaterials) {
        char fullname[4096];
        if (strlen(dirname))
            sprintf(fullname, "%s%s", dirname, mfilename);
        else
            sprintf(fullname, "%s", mfilename);
        if (!yo__save_mtl(fullname, scene)) return false;
    }

    FILE* file = fopen(filename, "wt");
    if (!file) return false;

    if (scene->nmaterials) {
        fprintf(file, "mtllib %s\n", mfilename);
    }

    yo__vert voffset = { 1, 1, 1, 1, 1 };

    // write cameras and environments if extensions are enabled
    if (ext) {
        for (int cid = 0; cid < scene->ncameras; cid++) {
            yo_camera* cam = &scene->cameras[cid];
            yo__fwrite_str(file, "o", ((cam->name) ? cam->name : ""));
            yo__fwrite_float3(file, "v", cam->from);
            yo__fwrite_float3(file, "v", cam->to);
            yo__fwrite_float3(file, "vn", cam->up);
            yo__fwrite_float3(file, "vn", cam->up);
            yo__fwrite_float2(file, "vt",
                              (float[2]){ cam->aperture, cam->aperture });
            yo__fwrite_float2(file, "vt",
                              (float[2]){ cam->width, cam->height });
            yo__fwrite_objverts(file, "c", 2, (int[2]){ 0, 1 }, voffset, 3,
                                (yo__vert){ 1, 1, 1, 0, 0, 0 });
            voffset.pos += 2;
            voffset.norm += 2;
            voffset.texcoord += 2;
        }
        for (int cid = 0; cid < scene->nenvs; cid++) {
            yo_env* env = &scene->envs[cid];
            yo__fwrite_str(file, "o", ((env->name) ? env->name : ""));
            yo__fwrite_str(file, "usemtl", (env->matname) ? env->matname : 0);
            yo__fwrite_float3(file, "v", env->from);
            yo__fwrite_float3(file, "v", env->to);
            yo__fwrite_float3(file, "vn", env->up);
            yo__fwrite_float3(file, "vn", env->up);
            yo__fwrite_float2(file, "vt", (float[2]){ 0, 0 });
            yo__fwrite_float2(file, "vt", (float[2]){ 0, 0 });
            yo__fwrite_objverts(file, "e", 2, (int[2]){ 0, 1 }, voffset, 3,
                                (yo__vert){ 1, 1, 1, 0, 0, 0 });
            voffset.pos += 2;
            voffset.norm += 2;
            voffset.texcoord += 2;
        }
    }

    // write all shape data
    for (int sid = 0; sid < scene->nshapes; sid++) {
        yo_shape* shape = &scene->shapes[sid];
        // shape header (name, material)
        yo__fwrite_str(file, "o", ((shape->name) ? shape->name : ""));
        yo__fwrite_str(file, "usemtl", (shape->matname) ? shape->matname : 0);
        if (ext && shape->xform) yo__fwrite_float16(file, "xf", shape->xform);

        // shape vertices
        yo__vert vto_write = { (shape->pos) ? 1 : 0,
                               (shape->texcoord) ? 1 : 0,
                               (shape->norm) ? 1 : 0,
                               ((ext && shape->color) ? 1 : 0),
                               ((ext && shape->radius) ? 1 : 0),
                               0 };
        int nto_write = 0;
        for (int i = 0; i < ((ext) ? 5 : 3); i++)
            nto_write = (&vto_write.pos)[i] ? i + 1 : nto_write;
        for (int j = 0; j < shape->nverts; j++) {
            yo__fwrite_float3(file, "v ", shape->pos + 3 * j);
            if (vto_write.norm)
                yo__fwrite_float3(file, "vn", shape->norm + 3 * j);
            if (vto_write.texcoord)
                yo__fwrite_float2(file, "vt", shape->texcoord + 2 * j);
            if (ext && vto_write.color)
                yo__fwrite_float3(file, "vc", shape->color + 3 * j);
            if (ext && vto_write.radius)
                yo__fwrite_float(file, "vr", shape->radius[j]);
        }

        // shape elements
        switch (shape->etype) {
            case yo_etype_point:
            case yo_etype_line:
            case yo_etype_triangle:
            case yo_etype_quad: {
                const char* labels[5] = { 0, "p", "l", "f", "f" };
                int esize = shape->etype;
                const char* label = labels[shape->etype];
                for (int j = 0; j < shape->nelems; j++) {
                    int* f = shape->elem + j * esize;
                    yo__fwrite_objverts(file, label, esize, f, voffset,
                                        nto_write, vto_write);
                }
            } break;
            case yo_etype_polyline:
            case yo_etype_polygon: {
                const char* label =
                    (shape->etype == yo_etype_polyline) ? "l" : "f";
                for (int j = 0, e = 0; j < shape->nelems; j++) {
                    int esize = shape->elem[e++];
                    int* f = shape->elem + e;
                    yo__fwrite_objverts(file, label, esize, f, voffset,
                                        nto_write, vto_write);
                }
            } break;
            default: { assert(false); } break;
        }
        for (int i = 0; i < 5; i++)
            (&voffset.pos)[i] += ((&vto_write.pos)[i]) ? shape->nverts : 0;
    }

    fclose(file);

    return true;
}

// -----------------------------------------------------------------------------
// BINARY DUMP LOADING
// -----------------------------------------------------------------------------

// magic code fron binary dump
#define yo__binmagic 0xaf45e782

// binary dump fixed integer array
static inline bool
yo__fread_binintn(FILE* file, int* v, int n) {
    if (fread(v, sizeof(int), n, file) != n) return false;
    return true;
}

// binary dump int array of variable length
static inline bool
yo__fread_binintarray(FILE* file, int** v) {
    int nn = 0;
    if (fread(&nn, sizeof(int), 1, file) != 1) return false;
    if (nn) {
        *v = (int*)calloc(nn, sizeof(int));
        if (fread(*v, sizeof(int), nn, file) != nn) return false;
    } else {
        *v = 0;
    }
    return true;
}

// binary dump fixed float array
static inline bool
yo__fread_binfloatn(FILE* file, float* v, int n) {
    if (fread(v, sizeof(float), n, file) != n) return false;
    return true;
}

// binary dump float array of varying size
static inline bool
yo__fread_binfloatarray(FILE* file, float** v) {
    int nn = 0;
    if (fread(&nn, sizeof(int), 1, file) != 1) return false;
    if (nn) {
        *v = (float*)calloc(nn, sizeof(float));
        if (fread(*v, sizeof(float), nn, file) != nn) return false;
    } else {
        *v = 0;
    }
    return true;
}

// binary dump strings
static inline bool
yo__fread_binstr(FILE* file, char** s) {
    int nn = 0;
    if (fread(&nn, sizeof(int), 1, file) != 1) return false;
    if (nn) {
        *s = (char*)calloc(nn, sizeof(char));
        if (fread(*s, 1, nn, file) != nn) return false;
    } else {
        *s = 0;
    }
    return true;
}

//
// binary dump OBJ (note that material data is dumped in the same file)
//
YO_API yo_scene*
yo_load_objbin(const char* filename, bool ext) {
    FILE* file = fopen(filename, "rb");
    if (!file) return 0;

    int magic = 0;
    yo__fread_binintn(file, &magic, 1);
    if (magic != yo__binmagic) return false;

    yo_scene* scene = (yo_scene*)calloc(1, sizeof(yo_scene));

    yo__fread_binintn(file, &scene->ncameras, 1);
    scene->cameras = (yo_camera*)calloc(scene->ncameras, sizeof(yo_camera));
    for (int i = 0; i < scene->ncameras; i++) {
        yo_camera* cam = scene->cameras + i;
        yo__fread_binstr(file, &cam->name);
        yo__fread_binfloatn(file, cam->from, 3);
        yo__fread_binfloatn(file, cam->to, 3);
        yo__fread_binfloatn(file, cam->up, 3);
        yo__fread_binfloatn(file, &cam->width, 2);
        yo__fread_binfloatn(file, &cam->height, 2);
        yo__fread_binfloatn(file, &cam->aperture, 1);
    }

    yo__fread_binintn(file, &scene->nenvs, 1);
    scene->envs = (yo_env*)calloc(scene->nenvs, sizeof(yo_env));
    for (int i = 0; i < scene->nenvs; i++) {
        yo_env* env = scene->envs + i;
        yo__fread_binstr(file, &env->name);
        yo__fread_binstr(file, &env->matname);
        yo__fread_binfloatn(file, env->from, 3);
        yo__fread_binfloatn(file, env->to, 3);
        yo__fread_binfloatn(file, env->up, 3);
    }

    if (!ext) {
        scene->ncameras = 0;
        free(scene->cameras);
        scene->cameras = 0;
        scene->nenvs = 0;
        free(scene->envs);
        scene->envs = 0;
    }

    yo__fread_binintn(file, &scene->nmaterials, 1);
    scene->materials =
        (yo_material*)calloc(scene->nmaterials, sizeof(yo_material));
    for (int i = 0; i < scene->nmaterials; i++) {
        yo_material* mat = scene->materials + i;
        yo__fread_binstr(file, &mat->name);
        yo__fread_binintn(file, &mat->illum, 1);
        yo__fread_binfloatn(file, mat->ke, 3);
        yo__fread_binfloatn(file, mat->ka, 3);
        yo__fread_binfloatn(file, mat->kd, 3);
        yo__fread_binfloatn(file, mat->ks, 3);
        yo__fread_binfloatn(file, mat->kr, 3);
        yo__fread_binfloatn(file, mat->kt, 3);
        yo__fread_binfloatn(file, &mat->ns, 1);
        yo__fread_binfloatn(file, &mat->ior, 1);
        yo__fread_binfloatn(file, &mat->op, 1);
        yo__fread_binstr(file, &mat->ke_txt);
        yo__fread_binstr(file, &mat->ka_txt);
        yo__fread_binstr(file, &mat->kd_txt);
        yo__fread_binstr(file, &mat->ks_txt);
        yo__fread_binstr(file, &mat->kr_txt);
        yo__fread_binstr(file, &mat->kt_txt);
        yo__fread_binstr(file, &mat->ns_txt);
        yo__fread_binstr(file, &mat->op_txt);
        yo__fread_binstr(file, &mat->ior_txt);
        yo__fread_binstr(file, &mat->bump_txt);
        yo__fread_binstr(file, &mat->disp_txt);
    }

    yo__fread_binintn(file, &scene->nshapes, 1);
    scene->shapes = (yo_shape*)calloc(scene->nshapes, sizeof(yo_shape));
    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = scene->shapes + i;
        yo__fread_binstr(file, &shape->name);
        yo__fread_binstr(file, &shape->groupname);
        yo__fread_binstr(file, &shape->matname);
        yo__fread_binintn(file, &shape->nelems, 1);
        yo__fread_binintarray(file, &shape->elem);
        yo__fread_binintn(file, &shape->etype, 1);
        yo__fread_binintn(file, &shape->nverts, 1);
        yo__fread_binfloatarray(file, &shape->pos);
        yo__fread_binfloatarray(file, &shape->norm);
        yo__fread_binfloatarray(file, &shape->texcoord);
        yo__fread_binfloatarray(file, &shape->color);
        yo__fread_binfloatarray(file, &shape->radius);
        if (ext) {
            if (shape->color) {
                free(shape->color);
                shape->color = 0;
            }
            if (shape->radius) {
                free(shape->radius);
                shape->radius = 0;
            }
        }
    }

    fclose(file);

    // fix ids
    yo__set_matids(scene);
    yo__set_textures(scene);

    return scene;
}

// -----------------------------------------------------------------------------
// BINARY DUMP SAVING
// -----------------------------------------------------------------------------

// load fixed int array
static inline bool
yo__fwrite_binintn(FILE* file, const int* v, int n) {
    if (fwrite(v, sizeof(int), n, file) != n) return false;
    return true;
}

// load variable size int array
static inline bool
yo__fwrite_binintarray(FILE* file, const int* v, int n) {
    if (!v) n = 0;
    if (fwrite(&n, sizeof(int), 1, file) != 1) return false;
    if (n) {
        if (fwrite(v, sizeof(int), n, file) != n) return false;
    }
    return true;
}

// load fixed float array
static inline bool
yo__fwrite_binfloatn(FILE* file, const float* v, int n) {
    if (fwrite(v, sizeof(float), n, file) != n) return false;
    return true;
}

// load variable size float array
static inline bool
yo__fwrite_binfloatarray(FILE* file, const float* v, int n) {
    if (!v) n = 0;
    if (fwrite(&n, sizeof(int), 1, file) != 1) return false;
    if (fwrite(v, sizeof(float), n, file) != n) return false;
    return true;
}

// load string
static inline bool
yo__fwrite_binstr(FILE* file, const char* s) {
    int l = (s) ? (int)strlen(s) + 1 : 0;
    if (fwrite(&l, sizeof(int), 1, file) != 1) return false;
    if (l) {
        if (fwrite(s, l, 1, file) != 1) return false;
    }
    return true;
}

//
// load binary obj dump
//
YO_API bool
yo_save_objbin(const char* filename, const yo_scene* scene, bool ext) {
    FILE* file = fopen(filename, "wb");
    if (!file) return false;

    int magic = yo__binmagic;
    yo__fwrite_binintn(file, &magic, 1);

    if (ext) {
        yo__fwrite_binintn(file, &scene->ncameras, 1);
        for (int i = 0; i < scene->ncameras; i++) {
            yo_camera* cam = scene->cameras + i;
            yo__fwrite_binstr(file, cam->name);
            yo__fwrite_binfloatn(file, cam->from, 3);
            yo__fwrite_binfloatn(file, cam->to, 3);
            yo__fwrite_binfloatn(file, cam->up, 3);
            yo__fwrite_binfloatn(file, &cam->width, 2);
            yo__fwrite_binfloatn(file, &cam->height, 2);
            yo__fwrite_binfloatn(file, &cam->aperture, 1);
        }
        yo__fwrite_binintn(file, &scene->nenvs, 1);
        for (int i = 0; i < scene->nenvs; i++) {
            yo_env* env = scene->envs + i;
            yo__fwrite_binstr(file, env->name);
            yo__fwrite_binstr(file, env->matname);
            yo__fwrite_binfloatn(file, env->from, 3);
            yo__fwrite_binfloatn(file, env->to, 3);
            yo__fwrite_binfloatn(file, env->up, 3);
        }
    } else {
        int zero = 0;
        yo__fwrite_binintn(file, &zero, 1);
        yo__fwrite_binintn(file, &zero, 1);
    }

    yo__fwrite_binintn(file, &scene->nmaterials, 1);
    for (int i = 0; i < scene->nmaterials; i++) {
        yo_material* mat = scene->materials + i;
        yo__fwrite_binstr(file, mat->name);
        yo__fwrite_binintn(file, &mat->illum, 1);
        yo__fwrite_binfloatn(file, mat->ke, 3);
        yo__fwrite_binfloatn(file, mat->ka, 3);
        yo__fwrite_binfloatn(file, mat->kd, 3);
        yo__fwrite_binfloatn(file, mat->ks, 3);
        yo__fwrite_binfloatn(file, mat->kr, 3);
        yo__fwrite_binfloatn(file, mat->kt, 3);
        yo__fwrite_binfloatn(file, &mat->ns, 1);
        yo__fwrite_binfloatn(file, &mat->ior, 1);
        yo__fwrite_binfloatn(file, &mat->op, 1);
        yo__fwrite_binstr(file, mat->ke_txt);
        yo__fwrite_binstr(file, mat->ka_txt);
        yo__fwrite_binstr(file, mat->kd_txt);
        yo__fwrite_binstr(file, mat->ks_txt);
        yo__fwrite_binstr(file, mat->kr_txt);
        yo__fwrite_binstr(file, mat->kt_txt);
        yo__fwrite_binstr(file, mat->ns_txt);
        yo__fwrite_binstr(file, mat->op_txt);
        yo__fwrite_binstr(file, mat->ior_txt);
        yo__fwrite_binstr(file, mat->bump_txt);
        yo__fwrite_binstr(file, mat->disp_txt);
    }

    yo__fwrite_binintn(file, &scene->nshapes, 1);
    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = scene->shapes + i;
        yo__fwrite_binstr(file, shape->name);
        yo__fwrite_binstr(file, shape->groupname);
        yo__fwrite_binstr(file, shape->matname);
        yo__fwrite_binintn(file, &shape->nelems, 1);
        if (shape->etype != yo_etype_polyline ||
            shape->etype != yo_etype_polygon) {
            yo__fwrite_binintarray(file, shape->elem,
                                   shape->nelems * shape->etype);
        } else {
            int n = 0;
            for (int i = 0; i < shape->etype; i++) {
                int en = shape->elem[i];
                n += en;
                i += en;
            }
            yo__fwrite_binintarray(file, shape->elem, n);
        }
        yo__fwrite_binintn(file, &shape->etype, 1);
        yo__fwrite_binintn(file, &shape->nverts, 1);
        yo__fwrite_binfloatarray(file, shape->pos, shape->nverts * 3);
        yo__fwrite_binfloatarray(file, shape->norm, shape->nverts * 3);
        yo__fwrite_binfloatarray(file, shape->texcoord, shape->nverts * 2);
        if (ext) {
            yo__fwrite_binfloatarray(file, shape->color, shape->nverts * 3);
            yo__fwrite_binfloatarray(file, shape->radius, shape->nverts);
        } else {
            yo__fwrite_binfloatarray(file, 0, shape->nverts * 3);
            yo__fwrite_binfloatarray(file, 0, shape->nverts);
        }
    }

    fclose(file);

    return true;
}

// -----------------------------------------------------------------------------
// TEXTURE HANDLING
// -----------------------------------------------------------------------------

//
// handles texture loading using stb_image.h
//
#ifndef YO_NOIMG

// stb_images causes a lot of warning and issues when including,
// try to reduce them using pragmas
#ifndef STBI_INCLUDE_STB_IMAGE_H

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_STATIC

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif

#include "../yocto/ext/stb_image.h"

#pragma GCC diagnostic pop

#endif

//
// load texture data
//
YO_API void
yo_load_textures(yo_scene* scene, const char* filename, int req_comp) {
    stbi_set_flip_vertically_on_load(1);
    char fullname[4096];
    for (int i = 0; i < scene->ntextures; i++) {
        yo__split_path(filename, fullname, 0, 0);
        strcat(fullname, scene->textures[i].path);
        scene->textures[i].pixels = stbi_loadf(
            fullname, &scene->textures[i].width, &scene->textures[i].height,
            &scene->textures[i].ncomp, req_comp);
    }
    stbi_set_flip_vertically_on_load(0);
}

#endif

#endif

#endif
