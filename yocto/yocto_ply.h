//
// # Yocto/Ply: Tiny library for PLY parsing/writing
//
// Yocto/Ply is a simple Stanford PLY parser that works with callbacks.
// The praser is designed for large files and does keep a copy of the model.
//
// Yocto/Ply also support writing PLY files again without keeping a copy of the
// model but instead writing elements directly after each call.
//
// Error reporting is done by throwing `std::runtime_error` exceptions.
//
// ## Parse a PLY file
//
// 1. define callbacks in `ply_callbacks` structure
// 2. run the parse with `load_ply()`
//
// ## Write amn PLY file
//
// 1. use `init_ply_stream()` to initialize the file streams for weriting
// 2. use the `write_ply_XXX()` function to write single Ply elements
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_PLY_H_
#define _YOCTO_PLY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <algorithm>

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Type of Ply data
enum struct ply_type { i8, i16, i32, i64, u8, u16, u32, u64, f32, f64 };

// Ply property
struct ply_property {
  string   name       = "";
  bool     is_list    = false;
  ply_type value_type = ply_type::f32;
  ply_type list_type  = ply_type::f32;
};

// Ply elements
struct ply_element {
  string               name       = "";
  int64_t              count      = 0;
  vector<ply_property> properties = {};
};

// Ply stream
struct ply_file {
  // Move-only object with automatic file closing on destruction.
  ply_file() {}
  ply_file(const string& filename) { open(filename); }
  ply_file(const ply_file&) = delete;
  ply_file& operator=(const ply_file&) = delete;
  ~ply_file() { close(); }

  // Open/Close file
  void open(const string& filename) {
    this->filename = filename;
    fs             = fopen(filename.c_str(), "rb");
    if (!fs) throw std::runtime_error{"cannot open file " + filename};
  }
  void close() { fclose(fs); }

  // Private data
  string filename   = "";
  FILE*  fs         = nullptr;
  bool   ascii      = false;
  bool   big_endian = false;
};

// Read Ply functions
void read_ply_header(ply_file& ply, vector<ply_element>& elements);
void read_ply_value(ply_file& ply, const ply_element& element,
    vector<double>& values, vector<vector<double>>& lists);
void read_ply_value(ply_file& ply, const ply_element& element,
    vector<float>& values, vector<vector<int>>& lists);

// Helpers to get element and property indices
int   find_ply_element(const vector<ply_element>& elements, const string& name);
int   find_ply_property(const ply_element& element, const string& name);
vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2);
vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3);
vec4i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3, const string& name4);

}  // namespace yocto

#endif
