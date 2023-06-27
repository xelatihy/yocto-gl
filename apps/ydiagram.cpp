//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include <yocto/yocto_cli.h>
#include <yocto/yocto_diagram.h>
#include <yocto/yocto_sceneio.h>

#include <functional>

namespace yocto {

using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto group       = "all"s;
  auto interactive = false;
  auto draft       = false;

  // diagrams and tags
  using diagram_func = std::function<diagram_data(const string&)>;
  auto diagrams = unordered_map<string, pair<vector<string>, diagram_func>>{
      // placeholder -------------------------------------------
      {"placeholder", {placeholder_tags, placeholder_diagrams}},
      // math --------------------------------------------------
      {"frame", {frame_tags, frame_diagrams}},
      {"objtransform", {objtransform_tags, objtransform_diagrams}},
      {"primitive", {primitive_tags, primitive_diagrams}},
      {"transform", {transform_tags, transform_diagrams}},
      {"vector", {vector_tags, vector_diagrams}},
      // image -------------------------------------------------
      {"compositing", {compositing_tags, compositing_diagrams}},
      {"image", {image_tags, image_diagrams}},
      {"tonemapping", {tonemapping_tags, tonemapping_diagrams}},
      {"tonemapplot", {tonemapplot_tags, tonemapplot_diagrams}},
      // scene -------------------------------------------------
      {"barycentric", {barycentric_tags, barycentric_diagrams}},
      {"camera", {camera_tags, camera_diagrams}},
      {"environment", {environment_tags, environment_diagrams}},
      {"shape", {shape_tags, shape_diagrams}},
      {"shapeapprox", {shapeapprox_tags, shapeapprox_diagrams}},
      {"texcoords", {texcoords_tags, texcoords_diagrams}},
      // raytracing --------------------------------------------
      {"antialiasing", {antialiasing_tags, antialiasing_diagrams}},
      {"brdfframe", {brdfframe_tags, brdfframe_diagrams}},
      {"brdfplot", {brdfplot_tags, brdfplot_diagrams}},
      {"cameraray", {cameraray_tags, cameraray_diagrams}},
      {"lighting", {lighting_tags, lighting_diagrams}},
      {"rendering", {rendering_tags, rendering_diagrams}},
      // intersect ----------------------------------------------
      {"bbox", {bbox_tags, bbox_diagrams}},
      {"intersect", {intersect_tags, intersect_diagrams}},
      // modeling -----------------------------------------------
      {"bezier", {bezier_tags, bezier_diagrams}},
      {"subcurve", {subcurve_tags, subcurve_diagrams}},
      {"subdiv", {subdiv_tags, subdiv_diagrams}},
      // renderingeq --------------------------------------------
      {"integration", {integration_tags, integration_diagrams}},
      {"mcplot", {mcplot_tags, mcplot_diagrams}},
      {"pisamples", {pisamples_tags, pisamples_diagrams}},
      {"sampling", {sampling_tags, sampling_diagrams}},
  };

  // labels
  auto labels = vector<string>{};
  for (auto& [group, _] : diagrams) labels.push_back(group);

  // parse command line
  auto cli = make_cli("ydiagram", "draw example diagrams");
  add_option(cli, "group", group, "Diagram group.");
  add_option(cli, "interactive", interactive, "Diagram interactive.");
  add_option(cli, "draft", draft, "Draft render.");
  parse_cli(cli, args);

  // render diagrams
  auto groups = group == "all" ? labels : vector<string>{group};
  for (auto& group : groups) {
    auto& [tags, diagrams_func] = diagrams.at(group);
    for (auto& tag : tags) {
      auto name = group + "_" + tag;
      print_info("rendering {}...", name);
      auto diagram = diagrams_func(tag);
      if (!interactive) {
        save_image("testsd/" + name + ".png", render_diagram(diagram, draft));
      } else {
        view_diagram(name, diagram);
      }
    }
  }
}

}  // namespace yocto

// Main
int main(int argc, const char* argv[]) {
  try {
    yocto::run({argv, argv + argc});
    return 0;
  } catch (const std::exception& error) {
    yocto::print_error(error.what());
    return 1;
  }
}
