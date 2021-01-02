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
// IMAGE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image view
struct imageview_state;
unique_ptr<imageview_state> make_imageview(const string& title);

// Run view
void run_view(imageview_state* viewer);

// Set image
void set_image(imageview_state* viewer, const string& name,
    const image<vec4f>& img, float exposure = 0, bool filmic = false);
void set_image(
    imageview_state* viewer, const string& name, const image<vec4b>& img);
void close_image(imageview_state* viewer, const string& name);

// Set params
void set_widget(imageview_state* viewer, const string& name,
    const string& pname, const json_value& param, const json_value& schema);
void set_widgets(imageview_state* viewer, const string& name,
    const json_value& params, const json_value& schema);

// Set ui callback
using imageview_callback =
    function<void(const string&, const json_value&, const gui_input&)>;
void set_callback(imageview_state* viewer, const imageview_callback& callback);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Make a scene view
struct traceview_state;
unique_ptr<traceview_state> make_traceview(const string& title);

// Run view
void run_view(traceview_state* viewer);

// Set scene and tracee params
void set_scene(
    imageview_state* viewer, const string& name, const trace_scene* scene);
void set_camera(
    imageview_state* viewer, const string& name, const trace_camera* camera);
void set_params(
    imageview_state* viewer, const string& name, const trace_params& params);
void close_scene(imageview_state* viewer, const string& name);

// Set ui callback
using traceview_callback =
    function<void(const string&, const json_value&, const gui_input&)>;
void set_callback(traceview_state* viewer, const traceview_callback& callback);

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
struct imageview_input {
  string name = "";

  bool close = false;

  bool         ichanged = true;
  image<vec4f> hdr      = {};
  image<vec4b> ldr      = {};
  float        exposure = 0;
  bool         filmic   = false;

  bool       wchanged = true;
  json_value widgets  = {};
  json_value schema   = {};
};

// An image visualized
struct imageview_view {
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

  // user params
  json_value widgets = json_value::object();
  json_value schema  = to_schema_object("User params.");

  ~imageview_view() {
    if (glimage) delete glimage;
  }
};

// Image pointer
using imageview_viewptr  = unique_ptr<imageview_view>;
using imageview_inputptr = unique_ptr<imageview_input>;

// Simple image viewer
struct imageview_state {
  vector<imageview_viewptr>  views       = {};       // views
  imageview_view*            selected    = nullptr;  // selected
  std::mutex                 input_mutex = {};
  vector<imageview_inputptr> inputs      = {};  // input images
  imageview_callback         callback    = {};  // params and ui callback
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE VIEWER
// -----------------------------------------------------------------------------
namespace yocto {

// Input image
struct traceview_input {
  string name = "";

  bool close = false;

  bool               schanged = true;
  const trace_scene* scene    = nullptr;

  bool                pchanged = true;
  const trace_camera* camera   = nullptr;
  trace_params        params   = {};

  bool       wchanged = true;
  json_value widgets  = {};
  json_value schema   = {};
};

// An image visualized
struct traceview_view {
  // original data
  string name = "scene.json";

  // trace data
  image<vec4f> hdr = {};
  image<vec4b> ldr = {};

  // diplay data
  image<vec4b> display  = {};
  float        exposure = 0;
  bool         filmic   = false;

  // viewing properties
  ogl_image*       glimage  = new ogl_image{};
  ogl_image_params glparams = {};

  // user params
  json_value params = json_value::object();
  json_value schema = to_schema_object("User params.");

  ~traceview_view() {
    if (glimage) delete glimage;
  }
};

// Image pointer
using traceview_viewptr  = unique_ptr<traceview_view>;
using traceview_inputptr = unique_ptr<traceview_input>;

// Simple image viewer
struct traceview_state {
  vector<traceview_viewptr>  images      = {};       // images
  traceview_view*            selected    = nullptr;  // selected
  std::mutex                 input_mutex = {};
  vector<traceview_inputptr> inputs      = {};  // input images
  imageview_callback         callback    = {};  // params and ui callback
};

}  // namespace yocto

#endif
