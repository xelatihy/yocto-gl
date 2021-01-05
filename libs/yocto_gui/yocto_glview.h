//
// Yocto/ImageViewer: Simpler image viewer.
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
//

#ifndef _YOCTO_IMAGEVIEWER_
#define _YOCTO_IMAGEVIEWER_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_json.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>
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
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void view_image(const image<vec4f>& img);
void view_image(const image<vec4b>& img);
void view_image(const image_data& img);

// Open a window and show a shape via path tracing
void view_shape(const string& title, const string& name,
    const shape_data& shape, bool addsky = false,
    const progress_callback& progress_cb = {});

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    const string& camera = "", const progress_callback& progress_cb = {});

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    scene_camera& camera, const trace_params& params,
    const progress_callback& progress_cb = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image view
struct ogl_imageviewer;
ogl_imageviewer make_imageviewer(const string& title);

// Run view
void run_viewer(ogl_imageviewer& viewer);

// Set image
void set_image(ogl_imageviewer& viewer, const string& name,
    const image<vec4f>& img, float exposure = 0, bool filmic = false);
void set_image(
    ogl_imageviewer& viewer, const string& name, const image<vec4b>& img);
void set_image(
    ogl_imageviewer& viewer, const string& name, const image_data& image);
void close_image(ogl_imageviewer& viewer, const string& name);

// Set params
void set_widget(ogl_imageviewer& viewer, const string& name,
    const string& pname, const json_value& param, const json_value& schema);
void set_widgets(ogl_imageviewer& viewer, const string& name,
    const json_value& params, const json_value& schema);

// Set ui callback
using ogl_imageviewer_callback =
    function<void(const string&, const json_value&, const gui_input&)>;
void set_callback(
    ogl_imageviewer& viewer, const ogl_imageviewer_callback& callback);

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

// Input image
struct ogl_imageinput {
  string name = "";

  bool close = false;

  bool       ichanged = true;
  image_data image    = {};
  float      exposure = 0;
  bool       filmic   = false;

  bool       wchanged = true;
  json_value widgets  = {};
  json_value schema   = {};
};

// An image visualized
struct ogl_imageview {
  // original data
  string name = "image.png";

  // image data
  image_data image = {};

  // diplay data
  image_data display  = {};
  float      exposure = 0;
  bool       filmic   = false;

  // viewing properties
  ogl_image        glimage  = {};
  ogl_image_params glparams = {};

  // user params
  json_value widgets = json_value::object();
  json_value schema  = to_schema_object("User params.");
};

// Image pointer
using ogl_imageview_ptr  = unique_ptr<ogl_imageview>;
using ogl_imageinput_ptr = unique_ptr<ogl_imageinput>;

// Simple image viewer
struct ogl_imageviewer {
  vector<ogl_imageview_ptr>  views       = {};       // views
  ogl_imageview*             selected    = nullptr;  // selected
  std::mutex                 input_mutex = {};
  vector<ogl_imageinput_ptr> inputs      = {};  // input images
  ogl_imageviewer_callback   callback    = {};  // params and ui callback
};

}  // namespace yocto

#endif
