//
// Yocto/ImageViewer: Simpler image viewer.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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
//

#ifndef _YOCTO_IMAGEVIEWER_
#define _YOCTO_IMAGEVIEWER_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_parallel.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::unique_ptr;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image view
struct imageview_state;
unique_ptr<imageview_state> make_imageview(const string& title);

// Run view
void run_view(imageview_state* viewer);

// Close viewer
void close_view(imageview_state* viewer);

// Set image
void set_image(imageview_state* viewer, const string& name,
    const image<vec4f>& img, float exposure = 0, bool filmic = false);
void set_image(
    imageview_state* viewer, const string& name, const image<vec4b>& img);
void close_image(imageview_state* viewer, const string& name);

// Open and asycn viewer
struct imageview_state;
unique_ptr<imageview_state> open_imageview(const string& title);
// Wait for the viewer to close
void wait_view(imageview_state* viewer);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Image viewer commands
enum imageview_command_type { quit, set, tonemap, close };
struct imageview_command {
  imageview_command_type type     = imageview_command_type::quit;
  string                 name     = "";
  image<vec4f>           hdr      = {};
  image<vec4b>           ldr      = {};
  float                  exposure = 0;
  bool                   filmic;
};

// Image view command queue and runner
using imageview_queue  = concurrent_queue<imageview_command>;
using imageview_runner = future<void>;

// An image visualized
struct imageview_image {
  // original data
  string name = "image.png";

  // image data
  image<vec4f> hdr = {};
  image<vec4b> ldr = {};

  // diplay data
  image<vec4b> display  = {};
  float        exposure = 0;
  bool         filmic   = false;

  // viewing properties
  ogl_image*       glimage  = new ogl_image{};
  ogl_image_params glparams = {};

  ~imageview_image() {
    if (glimage) delete glimage;
  }
};

// Image pointer
using imageview_imageptr = unique_ptr<imageview_image>;

// Simple image viewer
struct imageview_state {
  imageview_runner           runner;              // running thread
  imageview_queue            queue;               // command queue
  vector<imageview_imageptr> images   = {};       // images
  imageview_image*           selected = nullptr;  // selected
};

}  // namespace yocto

#endif
