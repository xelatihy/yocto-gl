//
// # Yocto/ShapeIO: Shape serialization
//
// Yocto/ShapeIO supports loading and saving shapes from Ply, Obj, Stl.
// Yocto/ShapeIO is implemented in `yocto_shapeio.h` and `yocto_shapeio.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#ifndef _YOCTO_SHAPEIO_H_
#define _YOCTO_SHAPEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
shape_data load_shape(const string& filename, bool flip_texcoords = true);
void       load_shape(
          const string& filename, shape_data& shape, bool flip_texcoords = true);
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
fvshape_data load_fvshape(const string& filename, bool flip_texcoords = true);
void         load_fvshape(
            const string& filename, fvshape_data& shape, bool flip_texcoords = true);
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
shape_data   make_shape_preset(const string& type);
fvshape_data make_fvshape_preset(const string& type);

// Load/save a shape
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
bool make_shape_preset(shape_data& shape, const string& type, string& error);
bool make_fvshape_preset(
    fvshape_data& shape, const string& type, string& error);

}  // namespace yocto

#endif
