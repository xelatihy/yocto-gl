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
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_shade.h>

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
void view_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show a set of images
void view_images(const string& title, const vector<string>& names,
    const vector<color_image>& images);

// Open a window and show an image for color grading
void colorgrade_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show a shape via path tracing
void view_shape(const string& title, const string& name,
    const scene_shape& shape, bool addsky = false);

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    bool print = true);

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    const string& camname, bool print = true);

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    const trace_params& params, bool print = true, bool edit = false);

using glview_scene_callback = std::function<void(gui_window* win,
    const gui_input& input, scene_model& scene, shade_scene& glscene)>;

void glview_scene(scene_model& scene, const string& name, const string& camname,
    const glview_scene_callback& widgets_callback  = {},
    const glview_scene_callback& uiupdate_callback = {});

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
void set_image(
    ogl_imageviewer& viewer, const string& name, const color_image& image);
void close_image(ogl_imageviewer& viewer, const string& name);

// Set params
void set_param(ogl_imageviewer& viewer, const string& name, const string& pname,
    const gui_param& param);
void set_params(ogl_imageviewer& viewer, const string& name,
    const string& pname, const gui_params& params);

// Set ui callback
using ogl_imageviewer_pcallback =
    function<void(const string&, const gui_params&)>;
using ogl_imageviewer_icallback =
    function<void(const string&, const gui_input&)>;
void set_params_callback(
    ogl_imageviewer& viewer, const ogl_imageviewer_pcallback& callback);
void set_input_callback(
    ogl_imageviewer& viewer, const ogl_imageviewer_icallback& callback);

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

  bool        ichanged = true;
  color_image image    = {};
  float       exposure = 0;
  bool        filmic   = false;

  bool       pchanged = true;
  string     pname    = "Params";
  gui_params params   = {};
};

// An image visualized
struct ogl_imageview {
  // original data
  string name = "image.png";

  // image data
  color_image image = {};

  // diplay data
  color_image display  = {};
  float       exposure = 0;
  bool        filmic   = false;

  // viewing properties
  ogl_image        glimage  = {};
  ogl_image_params glparams = {};

  // user params
  string     pname  = "Params";
  gui_params params = gui_params{};
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
  ogl_imageviewer_pcallback  pcallback   = {};  // params callback
  ogl_imageviewer_icallback  icallback   = {};  // input callback
};

}  // namespace yocto

#endif
