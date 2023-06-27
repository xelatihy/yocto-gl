//
// Implementation for Yocto/Diagram
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_diagram.h"

#include <yocto/yocto_sceneio.h>

#include <unordered_set>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM INITIALIZATION AND RENDERING
// -----------------------------------------------------------------------------
namespace yocto {

// Render text
texture_data make_text_texture(const string& text) {
  // Rendered text location
  static string text_location = "./textcache";

  // Latex template
  static string latex_template_header =
      "\\documentclass[10pt]{article}\n"
      "\\usepackage[cochineal]{newtx}\n"
      "\\usepackage[T1]{fontenc}\n"
      "\\usepackage{parskip}\n"
      "\n"
      "\\begin{document}\n"
      "\\pagenumbering{gobble}\n\n";
  static string latex_template_footer =
      "\n\n"
      "\\end{document}\n\n";

  // djb2 algorithm from http://www.cse.yorku.ca/~oz/hash.html
  auto hash = [](const string& text) -> string {
    auto str  = text.c_str();
    auto hash = (unsigned long)5381;
    auto c    = (int)*str++;
    while (c) {
      hash = ((hash << 5) + hash) + c;  // hash * 33 + c
      c    = (int)*str++;
    }
    return std::to_string(hash);
  };

  auto text_name    = hash(text);
  auto texture_name = text_location + "/" + text_name + ".png";
  auto latex_name   = text_location + "/" + text_name + ".tex";
  if (!path_exists(texture_name)) {
    // save latex
    printf("missing text: %s\n", text.c_str());
    auto latex = latex_template_header + text + latex_template_footer;
    save_text(latex_name, latex);
    // render latex
    auto command = "cd " + text_location + " && pdflatex " + text_name + ".tex";
    printf("%s\n", command.c_str());
    system(command.c_str());
    command = "cd " + text_location + " && pdfcrop " + text_name + ".pdf " +
              text_name + "-cropped.pdf";
    printf("%s\n", command.c_str());
    system(command.c_str());
    command = "cd " + text_location +
              " && gs -sDEVICE=pnggray -sBATCH -sOutputFile=" + text_name +
              ".png" + " -dNOPAUSE -r1200 " + text_name + "-cropped.pdf";
    printf("%s\n", command.c_str());
    system(command.c_str());
    system(("rm " + text_location + "/*.aux").c_str());
    system(("rm " + text_location + "/*.fls").c_str());
    system(("rm " + text_location + "/*.log").c_str());
    system(("rm " + text_location + "/*.gz").c_str());
    system(("rm " + text_location + "/*.fdb_latexmk").c_str());
  }

  auto texture = load_texture(texture_name);
  for (auto& c : texture.pixelsb) {
    c = {255, 255, 255, (byte)(255 - c.x)};
  }
  return texture;
}

// Text shape
shape_data make_text_shape(const string& text, const vec3f& offset,
    const texture_data& texture, float scale) {
  auto extents = max(texture.pixelsf.extents(), texture.pixelsb.extents());
  auto [width, height] = scale * extents / 72.0f;
  auto offset_         = vec3f{scale * 144 / 72.0f, scale * 144 / 72.0f, 1};
  auto z               = 0.1f;  // small offset to avoid z fighting
  auto shape           = shape_data{};
  shape.quads          = {{0, 1, 2, 3}};
  shape.positions      = {{-width, -height, z}, {+width, -height, z},
           {+width, +height, z}, {-width, +height, z}};
  for (auto& p : shape.positions) {
    if (offset.x > 0) p.x += +width + offset_.x * offset.x;
    if (offset.x < 0) p.x += -width + offset_.x * offset.x;
    if (offset.y > 0) p.y += +height + offset_.y * offset.y;
    if (offset.y < 0) p.y += -height + offset_.y * offset.y;
    if (offset.z > 0) p.z += offset.z;
  }
  shape.texcoords = {{0, 1}, {1, 1}, {1, 0}, {0, 0}};
  return shape;
}

// Text offset
frame3f make_text_frame(const string& text, const texture_data& texture,
    const vec3f& camera, const vec3f& position, float scale) {
  return lookat_frame(position, camera, {0, 1, 0}, true);
}

// Make text material
material_data make_text_material(const style_data& style, int textureid) {
  return {
      .type      = material_type::matte,
      .emission  = {style.text.x, style.text.y, style.text.z},
      .opacity   = style.text.w,
      .color_tex = textureid,
  };
}

// Make text positions
vector<vec3f> make_text_positions(
    const shape_data& shape, const label_data& text) {
  auto is_element = [](const string& label) {
    if (label.find("!!") == string::npos) return false;
    return label.substr(label.find("!!") + 2).find("e") != string::npos;
  };
  auto is_center = [](const string& label) {
    if (label.find("!!") == string::npos) return false;
    return label.substr(label.find("!!") + 2).find("c") != string::npos;
  };

  auto bbox = invalidb3f;
  for (auto& position : shape.positions) bbox = merge(bbox, position);
  auto center = (bbox.min + bbox.max) / 2;

  auto positions = vector<vec3f>{};
  for (auto&& [idx, label] : enumerate(text.labels)) {
    if (is_center(label)) {
      positions.push_back(center);
    } else if (!is_element(label)) {
      positions.push_back(text.positions.empty()
                              ? shape.positions[idx % shape.positions.size()]
                              : text.positions[idx % text.positions.size()]);
    } else {
      if (!shape.quads.empty()) {
        auto quad = shape.quads[idx % shape.quads.size()];
        positions.push_back(
            (shape.positions[quad.x] + shape.positions[quad.y] +
                shape.positions[quad.z] + shape.positions[quad.w]) /
            4);
      } else if (!shape.triangles.empty()) {
        auto triangle = shape.triangles[idx % shape.triangles.size()];
        positions.push_back(
            (shape.positions[triangle.x] + shape.positions[triangle.y] +
                shape.positions[triangle.z]) /
            3);
      } else if (!shape.lines.empty()) {
        auto line = shape.lines[idx % shape.lines.size()];
        positions.push_back(
            (shape.positions[line.x] + shape.positions[line.y]) / 2);
      } else if (!shape.points.empty()) {
        auto point = shape.points[idx % shape.points.size()];
        positions.push_back(shape.positions[point]);
      } else {
        positions.push_back({0, 0, 0});
      }
    }
  }
  return positions;
}

// split label
static pair<string, vec3f> split_label(const string& label) {
  auto pos = label.find("!!");
  if (pos == string::npos) return {label, {0, 0, 0}};
  auto text   = label.substr(0, pos);
  auto offset = vec3f{0, 0, 0};
  auto scale  = 1.0f;
  for (auto c : label.substr(pos + 2)) {
    auto reset_scale = true;
    switch (c) {
      case 'h':
        scale       = 0.5f;
        reset_scale = false;
        break;
      case 'a':
        scale       = 0.001f;
        reset_scale = false;
        break;
      case 'l': offset.x -= scale; break;
      case 'r': offset.x += scale; break;
      case 't': offset.y += scale; break;
      case 'b': offset.y -= scale; break;
      case 'n': offset.z += scale; break;
      case 'f': offset.z -= scale; break;
    }
    if (reset_scale) scale = 1.0f;
  }
  return {text, offset};
}

// Make text labels
vector<string> make_text_labels(const label_data& text) {
  if (text.labels.empty()) return {};
  auto labels = vector<string>{};
  for (auto& label : text.labels) labels.push_back(split_label(label).first);
  return labels;
}

// Make text offsets
vector<vec3f> make_text_offsets(const label_data& text) {
  if (text.labels.empty()) return {};
  auto offsets = vector<vec3f>{};
  for (auto& label : text.labels) offsets.push_back(split_label(label).second);
  return offsets;
}

// Make surface
shape_data make_fcolors_shape(
    const shape_data& shape, const vector<vec4f>& fcolors) {
  auto nshape = shape_data{};
  for (auto&& [fid, triangle] : enumerate(shape.triangles)) {
    auto offset = (int)nshape.positions.size();
    nshape.triangles.push_back({offset + 0, offset + 1, offset + 2});
    for (auto vid : triangle) {
      nshape.positions.push_back(shape.positions[vid]);
      if (!shape.texcoords.empty())
        nshape.texcoords.push_back(shape.texcoords[vid]);
      nshape.colors.push_back(fcolors[fid]);
    }
  }
  for (auto&& [fid, quad] : enumerate(shape.quads)) {
    auto offset = (int)nshape.positions.size();
    nshape.quads.push_back({offset + 0, offset + 1, offset + 2, offset + 3});
    for (auto vid : quad) {
      nshape.positions.push_back(shape.positions[vid]);
      if (!shape.texcoords.empty())
        nshape.texcoords.push_back(shape.texcoords[vid]);
      nshape.colors.push_back(fcolors[fid]);
    }
  }
  return nshape;
}

// Make surface
shape_data make_fill_shape(const shape_data& shape) {
  auto nshape      = shape_data{};
  nshape.triangles = shape.triangles;
  nshape.quads     = shape.quads;
  nshape.positions = shape.positions;
  nshape.texcoords = shape.texcoords;
  return nshape;
}

// Make fill material
material_data make_fill_material(const style_data& style, int textureid) {
  if (style.fcolors.empty()) {
    return {
        .type         = material_type::matte,
        .emission     = {style.fill.x, style.fill.y, style.fill.z},
        .opacity      = style.fill.w,
        .emission_tex = textureid,
    };
  } else {
    return {
        .type         = material_type::matte,
        .emission     = {1, 1, 1},
        .opacity      = 1,
        .emission_tex = textureid,
    };
  }
}

// Make a shape of lines and points
shape_data make_points_shape(
    const shape_data& shape, const frame3f& frame, float radius) {
  auto make_sphere = [](vec3f position, float radius) {
    auto frame = translation_frame(position) *
                 scaling_frame(vec3f{radius * 3, radius * 3, radius * 3});
    auto sphere = yocto::make_sphere(8);
    for (auto& p : sphere.positions) p = transform_point(frame, p);
    return sphere;
  };

  auto positions = shape.positions;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& point : shape.points) {
    merge_shape_inplace(nshape, make_sphere(positions[point], radius));
  }
  return nshape;
}

// Make a shape of lines
shape_data make_lines_shape(const shape_data& shape, const frame3f& frame,
    float radius, float connect) {
  auto make_capsule = [](vec3f position1, vec3f position2, float radius) {
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    return capsule;
  };

  auto positions = shape.positions;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& line : shape.lines) {
    if (connect <= 0) {
      merge_shape_inplace(
          nshape, make_capsule(positions[line.x], positions[line.y], radius));
    } else {
      auto position1 = positions[line.x], position2 = positions[line.y];
      auto middle = (position1 + position2) / 2;
      auto offset = max(0.0f, length(position2 - position1) - connect) *
                    normalize(position2 - position1) / 2;
      merge_shape_inplace(
          nshape, make_capsule(middle - offset, middle + offset, radius));
    }
  }
  return nshape;
}

// Make a shape of arrowa
shape_data make_arrows_shape(const shape_data& shape, const frame3f& frame,
    float radius, float connect) {
  auto make_capsule = [](vec3f position1, vec3f position2, float radius) {
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    return capsule;
  };
  auto make_arrow = [](vec3f position1, vec3f position2, float radius) {
    auto nshape         = shape_data{};
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    merge_shape_inplace(nshape, capsule);
    auto frame2 = translation_frame(position2) *
                  rotation_frame(rotation_axis, rotation_angle);
    auto cone = make_uvcone({8, 8, 8}, {radius * 3, radius * 5});
    for (auto& p : cone.positions) p = transform_point(frame2, p);
    merge_shape_inplace(nshape, cone);
    return nshape;
  };

  auto positions = shape.positions;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& line : shape.lines) {
    if (connect <= 0) {
      merge_shape_inplace(
          nshape, make_arrow(positions[line.x], positions[line.y], radius));
    } else {
      auto position1 = positions[line.x], position2 = positions[line.y];
      auto middle = (position1 + position2) / 2;
      auto offset = max(0.0f, length(position2 - position1) - connect) *
                    normalize(position2 - position1) / 2;
      merge_shape_inplace(
          nshape, make_arrow(middle - offset, middle + offset, radius));
    }
  }
  return nshape;
}

// Make stroke material
material_data make_stroke_material(
    const style_data& style, int textureid = invalidid) {
  return {
      .type         = material_type::matte,
      .emission     = {style.stroke.x, style.stroke.y, style.stroke.z},
      .opacity      = style.stroke.w,
      .emission_tex = textureid,
  };
}

// Make an image texture
texture_data make_image_texture(const array2d<vec4f>& image, bool nearest) {
  return texture_data{.pixelsb = float_to_byte(image), .nearest = nearest};
}

// Extract the wireframe
shape_data make_shape_wireframe(const shape_data& shape) {
  auto edges_triangles = get_edges(shape.triangles);
  auto edges_quads     = get_edges(shape.quads);
  auto edges           = vector<vec2i>{};
  edges.insert(edges.end(), edges_triangles.begin(), edges_triangles.end());
  edges.insert(edges.end(), edges_quads.begin(), edges_quads.end());
  return shape_data{
      .lines     = edges,
      .positions = shape.positions,
  };
}

// Render a diagram into a collection of scenes
void add_object(scene_data& scene, const frame3f& frame,
    const shape_data& shape, const label_data& text, const style_data& style) {
  // shape
  if ((!shape.triangles.empty() || !shape.quads.empty()) && style.fill.w > 0) {
    if (style.fcolors.empty()) {
      scene.shapes.push_back(make_fill_shape(shape));
    } else {
      scene.shapes.push_back(make_fcolors_shape(shape, style.fcolors));
    }
    if (style.texture.empty()) {
      scene.materials.push_back(make_fill_material(style, invalidid));
    } else {
      scene.textures.emplace_back(
          make_image_texture(style.texture, style.nearest));
      scene.materials.push_back(
          make_fill_material(style, (int)scene.textures.size() - 1));
    }
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = frame;
  }
  if (!shape.points.empty() && style.stroke.w > 0) {
    scene.shapes.push_back(make_points_shape(shape, frame, style.thickness));
    scene.materials.push_back(make_stroke_material(style));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }
  if (!shape.lines.empty() && !style.arrow && style.stroke.w > 0) {
    scene.shapes.push_back(
        make_lines_shape(shape, frame, style.thickness, style.connect));
    scene.materials.push_back(make_stroke_material(style));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }
  if (!shape.lines.empty() && style.arrow && style.stroke.w > 0) {
    scene.shapes.push_back(
        make_arrows_shape(shape, frame, style.thickness, style.connect));
    scene.materials.push_back(make_stroke_material(style));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }
  if (shape.lines.empty() &&
      (!shape.triangles.empty() || !shape.quads.empty()) && style.wireframe &&
      style.stroke.w > 0) {
    auto wireframe = make_shape_wireframe(shape);
    scene.shapes.push_back(
        make_lines_shape(wireframe, frame, style.thickness, style.connect));
    scene.materials.push_back(make_stroke_material(style));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // text
  if ((!text.labels.empty()) && style.text.w > 0) {
    auto positions = make_text_positions(shape, text);
    auto offsets   = make_text_offsets(text);
    auto labels    = make_text_labels(text);
    for (auto idx = 0; idx < (int)labels.size(); idx++) {
      if (labels[idx].empty()) continue;
      auto& texture = scene.textures.emplace_back(
          make_text_texture(labels[idx]));
      scene.materials.push_back(
          make_text_material(style, (int)scene.textures.size() - 1));
      scene.shapes.push_back(
          make_text_shape(labels[idx], offsets[idx], texture, style.textscale));
      auto& ninstance    = scene.instances.emplace_back();
      ninstance.frame    = make_text_frame(labels[idx], texture,
             scene.cameras.front().frame.o, transform_point(frame, positions[idx]),
             style.textscale);
      ninstance.material = (int)scene.materials.size() - 1;
      ninstance.shape    = (int)scene.shapes.size() - 1;
    }
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Make default camera
// 100 : 36 = 25 :  x => x = 36 * 25 / 100 = 9
// 100 : 36 =  d : 10 => d = 100 / 36 * 10 = 27.77777
//   l : 36 = 25 : 10 => l = 25 * 36 / 10 = 90
static camera_data default_camera() {
  return {
      .frame  = lookat_frame(vec3f{0, 0, 25}, vec3f{0, 0, 0}, vec3f{0, 1, 0}),
      .lens   = 0.09f,
      .aspect = 16.0f / 9.0f,
      .focus  = 25,
  };
}

// Init diagram
diagram_data make_diagram() {
  auto diagram = diagram_data{};
  diagram.scene.cameras.push_back(default_camera());
  return diagram;
}

// Group nodes
frame3f add_group(
    diagram_data& diagram, const vec3f& offset, const frame3f& frame) {
  return dtranslation(offset) * frame;
}
frame3f add_tgroup(diagram_data& diagram, const frame3f& frame) {
  return frame;
}
frame3f add_tgroup(
    diagram_data& diagram, const vec3f& offset, const frame3f& frame) {
  return dtranslation(offset) * frame;
}

// Shape nodes
void add_label(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, const string& label, const vec4f& textcolor) {
  add_object(diagram.scene, frame, {},
      {.positions = {position}, .labels = {label}}, {.text = textcolor});
}
void add_labels(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& textcolor) {
  add_object(diagram.scene, frame, {},
      {.positions = positions, .labels = labels}, {.text = textcolor});
}
void add_point(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, const string& label, const vec4f& stroke,
    float thickness) {
  add_points(diagram, frame, {position}, {label}, stroke, thickness);
}
void add_points(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto points = vector<int>{};
  for (auto idx : range((int)positions.size())) points.push_back(idx);
  add_object(diagram.scene, frame, {.points = points, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_points(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& positions_, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto positions = vector<vec3f>{};
  for (auto position : positions_)
    positions.push_back({position.x, position.y, 0});
  auto points = vector<int>{};
  for (auto idx : range((int)positions.size())) points.push_back(idx);
  add_object(diagram.scene, frame, {.points = points, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_line(diagram_data& diagram, const frame3f& frame,
    const vec3f& position1, const vec3f& position2,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  add_lines(diagram, frame, {position1, position2}, labels, stroke, thickness);
}
void add_lines(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto lines = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_lines(diagram_data& diagram, const frame3f& frame,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_arrow(diagram_data& diagram, const frame3f& frame,
    const vec3f& position1, const vec3f& position2,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  add_arrows(diagram, frame, {position1, position2}, labels, stroke, thickness);
}
void add_arrows(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto lines = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness, .arrow = true});
}
void add_vector(diagram_data& diagram, const frame3f& frame,
    const vec3f& direction, const vector<string>& labels, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame,
      {.lines = {{0, 1}}, .positions = {{0.0, 0.0, 0.0}, direction}},
      {.positions = {{0.0, 0.0, 0.0}, direction}, .labels = labels},
      {.stroke = stroke, .thickness = thickness, .arrow = true});
}
void add_axes(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  return add_axes(diagram, frame, translation_frame(position), scale, labels,
      stroke, thickness);
}
void add_axes(diagram_data& diagram, const frame3f& frame,
    const frame3f& frame1, float scale, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.lines        = {{0, 1}, {0, 2}, {0, 3}},
          .positions = {frame1.o, frame1.o + scale * frame1.x,
              frame1.o + scale * frame1.y, frame1.o + scale * frame1.z}},
      {.positions = {frame1.o, frame1.o + scale * frame1.x,
           frame1.o + scale * frame1.y, frame1.o + scale * frame1.z},
          .labels = labels},
      {.stroke = stroke, .thickness = thickness, .arrow = true});
}
void add_ray(diagram_data& diagram, const frame3f& frame, const vec3f& position,
    const vec3f& direction, const vector<string>& labels, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame,
      {.lines = {{0, 1}}, .positions = {position, position + direction}},
      {.positions = {position, position + direction}, .labels = labels},
      {.stroke = stroke, .thickness = thickness, .arrow = true});
}
void add_polyline(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto lines = vector<vec2i>{};
  for (auto idx : range((int)positions.size() - 1))
    lines.push_back({idx, idx + 1});
  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.positions = positions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_quad(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto positions = dconstants::quad_positions;
  for (auto& pos : positions) pos = pos * scale + position;
  add_object(diagram.scene, frame,
      {.quads = {{0, 1, 2, 3}}, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_quadv(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.quads = {{0, 1, 2, 3}}, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_tquad(diagram_data& diagram, const frame3f& frame,
    const array2d<vec4f>& texture, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.quads        = {{0, 1, 2, 3}},
          .positions = dconstants::quad_positions,
          .texcoords = dconstants::quad_texcoords},
      {},
      {.stroke       = stroke,
          .fill      = {1, 1, 1, 1},
          .thickness = thickness,
          .texture   = texture});
}
void add_rect(diagram_data& diagram, const frame3f& frame, const vec2f& aspect,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto positions = dconstants::quad_positions;
  for (auto& pos : positions)
    pos = pos * vec3f{aspect.x / aspect.y, 1, 1} * scale + position;
  add_object(diagram.scene, frame,
      {.quads = {{0, 1, 2, 3}}, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_triangle(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto positions = dconstants::triangle_positions;
  for (auto& pos : positions) pos = pos * scale + position;
  add_object(diagram.scene, frame,
      {.triangles = {{0, 1, 2}}, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_trianglev(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.triangles = {{0, 1, 2}}, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_polygon(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto center = vec3f{0, 0, 0};
  for (auto& position : positions) center += position / (float)positions.size();
  auto positions_ = vector<vec3f>{};
  positions_.push_back(center);
  positions_.insert(positions_.end(), positions.begin(), positions.end());
  auto triangles = vector<vec3i>{};
  auto lines     = vector<vec2i>{};
  for (auto idx : range((int)positions.size())) {
    triangles.push_back({0, idx + 1, (idx + 1) % (int)positions.size() + 1});
    lines.push_back({idx + 1, (idx + 1) % (int)positions.size() + 1});
  }
  add_object(diagram.scene, frame,
      {.lines = lines, .triangles = triangles, .positions = positions_},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_triangles(diagram_data& diagram, const frame3f& frame,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<string>& labels, const vec4f& fill, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame,
      {.triangles = triangles, .positions = positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_quads(diagram_data& diagram, const frame3f& frame,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<string>& labels, const vec4f& fill, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame, {.quads = quads, .positions = positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_disk(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto steps     = 32;
  auto positions = vector<vec3f>{};
  auto lines     = vector<vec2i>{};
  auto triangles = vector<vec3i>{};
  for (auto idx : range(steps)) {
    auto theta = 2 * pif * idx / (float)steps;
    positions.push_back(position + scale * vec3f{cos(theta), sin(theta), 0});
    triangles.push_back({idx, (idx + 1) % steps, steps});
    lines.push_back({idx, (idx + 1) % steps});
  }
  positions.push_back(position);
  add_object(diagram.scene, frame,
      {.lines = lines, .triangles = triangles, .positions = positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_circle(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto steps     = 32;
  auto positions = vector<vec3f>{};
  auto lines     = vector<vec2i>{};
  auto triangles = vector<vec3i>{};
  for (auto idx : range(steps)) {
    auto theta = 2 * pif * idx / (float)steps;
    positions.push_back(position + scale * vec3f{cos(theta), sin(theta), 0});
    triangles.push_back({idx, (idx + 1) % steps, steps});
    lines.push_back({idx, (idx + 1) % steps});
  }
  positions.push_back(position);
  add_object(diagram.scene, frame,
      {.lines = lines, .triangles = triangles, .positions = positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_arc(diagram_data& diagram, const frame3f& frame, const vec3f& position,
    float scale, float angle, const vector<string>& labels, const vec4f& stroke,
    float thickness) {
  auto steps     = 8;
  auto positions = vector<vec3f>{};
  for (auto idx : range(steps + 1)) {
    auto theta = angle * idx / (float)steps;
    positions.push_back(position + scale * vec3f{cos(theta), sin(theta), 0});
  }
  auto lines = vector<vec2i>{};
  for (auto idx : range(steps)) lines.push_back({idx, idx + 1});
  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_cube(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto shape = transform_shape(
      make_cube(), translation_frame(position) * scaling_frame(scale));
  add_object(diagram.scene, frame,
      {.quads = shape.quads, .positions = shape.positions}, {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_sphere(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto sphere = transform_shape(
      make_sphere(), translation_frame(position) * scaling_frame(scale));
  add_object(diagram.scene, frame,
      {.quads = sphere.quads, .positions = sphere.positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_uvsphere(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto sphere = transform_shape(
      make_uvsphere(), translation_frame(position) * scaling_frame(scale));
  add_object(diagram.scene, frame,
      {.quads = sphere.quads, .positions = sphere.positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_uvspheret(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& texture,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto sphere = transform_shape(
      make_sphere(), translation_frame(position) * scaling_frame(scale));
  add_object(diagram.scene, frame,
      {.quads = sphere.quads, .positions = sphere.positions},
      {.labels = labels},
      {.stroke       = stroke,
          .fill      = {1, 1, 1, 1},
          .thickness = thickness,
          .texture   = texture});
}
void add_cline(diagram_data& diagram, const frame3f& frame,
    const vec3f& position1, const vec3f& position2,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.lines = {{0, 1}}, .positions = {position1, position2}},
      {.positions = {position1, position2}, .labels = labels},
      {.stroke = stroke, .thickness = thickness, .connect = 0.3});
}
void add_carrow(diagram_data& diagram, const frame3f& frame,
    const vec3f& position1, const vec3f& position2,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame,
      {.lines = {{0, 1}}, .positions = {position1, position2}},
      {.positions = {position1, position2}, .labels = labels},
      {.stroke       = stroke,
          .thickness = thickness,
          .arrow     = true,
          .connect   = 0.3});
}
void add_grid(diagram_data& diagram, const frame3f& frame, const vec2i& steps,
    const vec4f& stroke, float thickness) {
  auto ratio     = (float)steps.x / (float)steps.y;
  auto lines     = vector<vec2i>{};
  auto positions = vector<vec3f>{};
  for (auto i = 0; i <= steps.x; i++) {
    auto u = (float)i / (float)steps.x;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back({(2 * u - 1) * ratio, -1, 0});
    positions.push_back({(2 * u - 1) * ratio, +1, 0});
  }
  for (auto i = 0; i <= steps.y; i++) {
    auto u = (float)i / (float)steps.y;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back({-ratio, 2 * u - 1, 0});
    positions.push_back({+ratio, 2 * u - 1, 0});
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions}, {},
      {.stroke = stroke, .thickness = thickness});
}
void add_grid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto ratio     = (float)steps.x / (float)steps.y;
  auto lines     = vector<vec2i>{};
  auto positions = vector<vec3f>{};
  for (auto i = 0; i <= steps.x; i++) {
    auto u = (float)i / (float)steps.x;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(position + scale * vec3f{(2 * u - 1) * ratio, -1, 0});
    positions.push_back(position + scale * vec3f{(2 * u - 1) * ratio, +1, 0});
  }
  for (auto i = 0; i <= steps.y; i++) {
    auto u = (float)i / (float)steps.y;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(position + scale * vec3f{-ratio, 2 * u - 1, 0});
    positions.push_back(position + scale * vec3f{+ratio, 2 * u - 1, 0});
  }
  auto lpositions = vector<vec3f>{};
  for (auto j = 0; j < steps.y; j++) {
    for (auto i = 0; i < steps.x; i++) {
      auto u = (i + 0.5f) / (float)steps.x;
      auto v = (j + 0.5f) / (float)steps.y;
      lpositions.push_back(
          position + scale * vec3f{(2 * u - 1) * ratio, 2 * v - 1, 0});
    }
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.positions = lpositions, .labels = labels},
      {.stroke = stroke, .thickness = thickness});
}
void add_cgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto csteps    = 32;
  auto lines     = vector<vec2i>{};
  auto positions = vector<vec3f>{};
  for (auto j = 1; j <= steps.x; j++) {
    auto radius = ((float)j / steps.x);
    for (auto i = 0; i < csteps; i++) {
      auto theta1 = 2 * pif * (float)(i + 0) / 32;
      auto theta2 = 2 * pif * (float)(i + 1) / 32;
      lines.push_back({(int)positions.size(), (int)positions.size() + 1});
      positions.push_back(
          position + scale * radius * vec3f{cos(theta1), sin(theta1), 0});
      positions.push_back(
          position + scale * radius * vec3f{cos(theta2), sin(theta2), 0});
    }
  }
  for (auto i = 0; i < steps.y; i++) {
    auto theta = 2 * pif * (float)i / steps.y;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(position);
    positions.push_back(position + scale * vec3f{cos(theta), sin(theta), 0});
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_cgridnu(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto csteps    = 32;
  auto lines     = vector<vec2i>{};
  auto positions = vector<vec3f>{};
  for (auto j = 1; j <= steps.x; j++) {
    auto radius = sqrt((float)j / steps.x);
    for (auto i = 0; i < csteps; i++) {
      auto theta1 = 2 * pif * (float)(i + 0) / 32;
      auto theta2 = 2 * pif * (float)(i + 1) / 32;
      lines.push_back({(int)positions.size(), (int)positions.size() + 1});
      positions.push_back(
          position + scale * radius * vec3f{cos(theta1), sin(theta1), 0});
      positions.push_back(
          position + scale * radius * vec3f{cos(theta2), sin(theta2), 0});
    }
  }
  for (auto i = 0; i < steps.y; i++) {
    auto theta = 2 * pif * (float)i / steps.y;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(position);
    positions.push_back(position + scale * vec3f{cos(theta), sin(theta), 0});
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_tgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto positions = dconstants::triangle_positions;
  for (auto& pos : positions) pos = position + pos * scale;
  return add_tgridv(
      diagram, frame, positions, steps, labels, stroke, thickness);
}
void add_tgridv(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& triangle, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto triangle_point = [&](const vec2f& uv) {
    return interpolate_triangle(
        triangle.at(0), triangle.at(1), triangle.at(2), sample_triangle(uv));
  };

  auto lines     = vector<vec2i>{};
  auto positions = vector<vec3f>{};
  for (auto i = 0; i <= steps.x; i++) {
    auto u = (float)i / steps.x;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(triangle_point({u, 0}));
    positions.push_back(triangle_point({u, 1}));
  }
  for (auto i = 0; i <= steps.y; i++) {
    auto v = (float)i / steps.y;
    lines.push_back({(int)positions.size(), (int)positions.size() + 1});
    positions.push_back(triangle_point({0, v}));
    positions.push_back(triangle_point({1, v}));
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_bezier(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& controlp, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto steps        = 64;
  auto bezier_point = [](vec3f p1, vec3f p2, vec3f p3, vec3f p4, float u) {
    return (1 - u) * (1 - u) * (1 - u) * p1 + 3 * (1 - u) * (1 - u) * u * p2 +
           3 * (1 - u) * u * u * p3 + u * u * u * p4;
  };
  auto positions = vector<vec3f>(steps + 1);
  for (auto idx : range(steps + 1)) {
    positions[idx] = bezier_point(
        controlp[0], controlp[1], controlp[2], controlp[3], idx / (float)steps);
  }
  return add_polyline(diagram, frame, positions, labels, stroke, thickness);
}
void add_qbezier(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& controlp, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto steps         = 64;
  auto qbezier_point = [](vec3f p1, vec3f p2, vec3f p3, float u) {
    return (1 - u) * (1 - u) * p1 + 2 * (1 - u) * u * p2 + u * u * p3;
  };
  auto positions = vector<vec3f>(steps + 1);
  for (auto idx : range(steps + 1)) {
    positions[idx] = qbezier_point(
        controlp[0], controlp[1], controlp[2], idx / (float)steps);
  }
  return add_polyline(diagram, frame, positions, labels, stroke, thickness);
}
void add_slabs(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto positions = dconstants::slabs_positions;
  for (auto& pos : positions) pos = position + pos * scale;
  add_object(diagram.scene, frame,
      {.quads = dconstants::slabs_quads, .positions = positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_slab(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels,
    const vec4f& fill, const vec4f& stroke, float thickness) {
  auto positions = dconstants::slab_positions;
  for (auto& pos : positions) pos = position + pos * scale;
  add_object(diagram.scene, frame,
      {.quads = dconstants::slab_quads, .positions = positions},
      {.labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_image(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& image,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto ratio     = (float)image.extent(0) / (float)image.extent(1);
  auto positions = dconstants::quad_positions;
  for (auto& pos : positions) pos = position + pos * scale * vec3f{ratio, 1, 1};
  add_object(diagram.scene, frame,
      {.quads        = dconstants::quad_quads,
          .positions = positions,
          .texcoords = dconstants::quad_texcoords},
      {.labels = labels},
      {.stroke       = stroke,
          .fill      = {1, 1, 1, 1},
          .thickness = thickness,
          .texture   = image,
          .nearest   = true});
}
void add_imagegrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& image,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto ratio = (float)image.extent(0) / (float)image.extent(1);
  add_image(diagram, frame, position, scale, image, {}, dcolors::transparent);
  add_grid(diagram, frame, position, scale, (vec2i)image.extents(), labels,
      stroke, thickness);
}
void add_imagefull(diagram_data& diagram, const frame3f& frame,
    const array2d<vec4f>& image, const string& label) {
  auto ratio = (float)image.extent(0) / (float)image.extent(1);
  add_tquad(diagram, dscaling({1, 1 / ratio, 1}) * dscaling(5), image,
      dcolors::transparent);
  if (!label.empty()) {
    add_rect(diagram, frame, {3, 1}, {3.8, -2.0, 0.1}, 0.24, {},
        dcolors::transparent, dcolors::black);
    add_label(diagram, frame, {3.8, -2.0, 0.1}, label, dcolors::white);
  }
}
void add_imagemosaic(diagram_data& diagram, const frame3f& frame,
    const vector<array2d<vec4f>>& images, const vector<string>& labels) {
  if (images.size() != 4)
    throw std::out_of_range{"supports only 4 images for now"};
  auto positions = vector<vec3f>{
      {-2.5, 1.25, 0}, {+2.5, 1.25, 0}, {-2.5, -1.25, 0}, {+2.5, -1.25, 0}};
  for (auto idx : range((int)size(images))) {
    auto& image     = images[idx];
    auto  transform = frame * dtranslation(positions[idx]);
    auto  ratio     = (float)image.extent(0) / (float)image.extent(1);
    add_tquad(diagram, transform * dscaling({1, 1 / ratio, 1}) * dscaling(2.5),
        image, dcolors::transparent);
    if (idx < labels.size() && !labels[idx].empty()) {
      add_rect(diagram, transform, {2, 1}, {1.9, -0.94, 0.1}, 0.175, {},
          dcolors::transparent, dcolors::black);
      add_label(
          diagram, transform, {1.9, -0.94, 0.1}, labels[idx], dcolors::white);
    }
  }
}
void add_rpoints(diagram_data& diagram, const frame3f& frame,
    const diagram_func<vec3f, vec2f>& func, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto rng       = make_rng(187291);
  auto positions = vector<vec3f>{};
  auto points    = vector<int>{};
  auto ssteps    = int(round(sqrt((float)steps)));
  for (auto idx : range(steps)) {
    auto uv = !stratified ? rand2f(rng)
                          : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                (idx / ssteps + rand1f(rng)) / ssteps};
    positions.push_back(func(uv));
    points.push_back(idx);
  }
  add_object(diagram.scene, frame, {.points = points, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_rpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_rpoints(
      diagram, frame,
      [&](vec2f uv) {
        return position + scale * vec3f{2 * uv.x - 1, 2 * uv.y - 1, 0};
      },
      steps, stratified, labels, stroke, thickness);
}
void add_rdpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_rpoints(
      diagram, frame,
      [&](vec2f uv) {
        return position +
               scale * vec3f{sample_disk(uv).x, sample_disk(uv).y, 0};
      },
      steps, stratified, labels, stroke, thickness);
}
void add_rdpointsnu(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_rpoints(
      diagram, frame,
      [&](vec2f uv_) {
        auto uv = vec2f{uv_.x, uv_.y * uv_.y};
        return position +
               scale * vec3f{sample_disk(uv).x, sample_disk(uv).y, 0};
      },
      steps, stratified, labels, stroke, thickness);
}
void add_rtpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto positions = dconstants::triangle_positions;
  for (auto& pos : positions) pos = position + pos * scale;
  return add_rtpoints(
      diagram, frame, positions, steps, stratified, labels, stroke, thickness);
}
void add_rtpoints(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& triangle, int steps, bool stratified,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_rpoints(
      diagram, frame,
      [&](vec2f uv_) {
        auto uv = sample_triangle(uv_);
        return interpolate_triangle(triangle[0], triangle[1], triangle[2], uv);
      },
      steps, stratified, labels, stroke, thickness);
}
void add_rlines(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, rlines_type type, int steps,
    bool stratified, const vector<string>& labels, const vec4f& stroke,
    float thickness) {
  auto rng       = make_rng(187291);
  auto ssteps    = int(round(sqrt((float)steps)));
  auto positions = vector<vec3f>{};
  auto lines     = vector<vec2i>{};
  for (auto idx : range(steps)) {
    auto uv = !stratified ? rand2f(rng)
                          : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                (idx / ssteps + rand1f(rng)) / ssteps};
    switch (type) {
      case rlines_type::hemi:
        positions.push_back(position);
        positions.push_back(
            position + scale * sample_hemisphere({0, 0, 1}, uv));
        break;
      case rlines_type::hemicos:
        positions.push_back(position);
        positions.push_back(
            position + scale * sample_hemisphere_cos({0, 0, 1}, uv));
        break;
      case rlines_type::hemicospower:
        positions.push_back(position);
        positions.push_back(
            position + scale * sample_hemisphere_cospower(64, {0, 0, 1}, uv));
        break;
      case rlines_type::beam:
        auto rdisk = sample_disk(uv);
        positions.push_back(position + scale * vec3f{rdisk.x, rdisk.y, 0});
        positions.push_back(position + scale * vec3f{rdisk.x, rdisk.y, 1});
        break;
    }
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  }

  add_object(diagram.scene, frame, {.lines = lines, .positions = positions},
      {.labels = labels}, {.stroke = stroke, .thickness = thickness});
}
void add_atgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  auto shape = shape_data{};
  for (auto i = -steps.x / 2; i <= steps.x; i++) {
    auto u = i / (float)steps.x;
    auto a = vec3f{-1 + u * 2 - 0.075f, -1.1f, 0},
         b = vec3f{u * 2 + 0.075f, 1.1f, 0};
    if (a.x < -1.1f) {
      auto t = (-1.1f - a.x) / (b.x - a.x);
      a      = a + t * (b - a);
    }
    if (b.x > +1.1f) {
      auto t = (+1.1f - a.x) / (b.x - a.x);
      b      = a + t * (b - a);
    }
    shape.lines.push_back(
        {(int)shape.positions.size(), (int)shape.positions.size() + 1});
    shape.positions.push_back(position + scale * a);
    shape.positions.push_back(position + scale * b);
  }
  for (auto i = 0; i <= steps.y; i++) {
    auto u = i / (float)steps.y;
    shape.lines.push_back(
        {(int)shape.positions.size(), (int)shape.positions.size() + 1});
    shape.positions.push_back(position + scale * vec3f{-1.1f, -1 + u * 2, 0});
    shape.positions.push_back(position + scale * vec3f{+1.1f, -1 + u * 2, 0});
  }

  add_object(diagram.scene, frame, shape, {.labels = labels},
      {.stroke = stroke, .thickness = thickness});
}

// Shape nodes
void add_lines(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame, shape, {},
      {.stroke = stroke, .fill = dcolors::transparent, .thickness = thickness});
}
void add_lines(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame, shape,
      {.positions = shape.positions, .labels = labels},
      {.stroke = stroke, .fill = dcolors::transparent, .thickness = thickness});
}
void add_shape(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vec4f& fill, const vec4f& stroke,
    float thickness) {
  add_object(diagram.scene, frame, shape, {},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_shape(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels, const vec4f& fill,
    const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame, shape,
      {.positions = shape.positions, .labels = labels},
      {.stroke = stroke, .fill = fill, .thickness = thickness});
}
void add_shape(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels,
    const array2d<vec4f>& texture, const vec4f& stroke, float thickness) {
  add_object(diagram.scene, frame, shape, {},
      {.stroke       = stroke,
          .fill      = {1, 1, 1, 1},
          .thickness = thickness,
          .texture   = texture});
}

// Easy checking in lists
bool contains(const string& tag, const vector<string>& tags) {
  for (auto tag_ : tags) {
    if (tag_ == tag) return true;
  }
  return false;
}

// Frames
frame3f dtranslation(const vec3f& translation) {
  return translation_frame(translation);
}
frame3f drotation(const vec3f& rotation) {
  return rotation_frame(vec3f{1, 0, 0}, radians(rotation.x)) *
         rotation_frame(vec3f{0, 1, 0}, radians(rotation.y)) *
         rotation_frame(vec3f{0, 0, 1}, radians(rotation.z));
}
frame3f dscaling(const vec3f& scaling) { return scaling_frame(scaling); }
frame3f dscaling(float scale) {
  return scaling_frame(vec3f{scale, scale, scale});
}
frame3f dlookat(const vec3f& from, const vec3f& to, const vec3f& up) {
  return lookat_frame(from, to, up, true);
}
frame3f dtransform(const vec3f& translation, const vec3f& rotation) {
  return dtranslation(translation) * drotation(rotation);
}
frame3f dgtransform(
    const vec3f& translation, const vec3f& rotation, const vec3f& scaling) {
  return dtranslation(translation) * drotation(rotation) * dscaling(scaling);
}

// Shape transformations
shape_data scale_shape(const shape_data& shape, float size) {
  auto tshape = shape;
  for (auto& position : tshape.positions) position *= size;
  return tshape;
}
shape_data transform_shape(
    const shape_data& shape, const frame3f& frame, float size) {
  return transform_shape(shape, frame * scaling_frame(vec3f{size, size, size}));
}

// Shapes
shape_data load_dshape(const string& filename) {
  auto loaded = load_shape(filename);
  auto shape  = shape_data{
       .points    = loaded.points,
       .lines     = loaded.lines,
       .triangles = loaded.triangles,
       .quads     = loaded.quads,
       .positions = loaded.positions,
       .texcoords = loaded.texcoords,
  };
  return shape;
}

// Plot functions
frame3f add_plot(
    diagram_data& diagram, const frame3f& frame, const dplotaxes_params& axes) {
  // add frame
  auto ratio = axes.aspect.x / axes.aspect.y;
  add_rect(diagram, frame, {ratio, 1}, {0, 0, 0}, 1, {}, dcolors::transparent,
      dcolors::black);
  for (auto& title : axes.titles) {
    add_label(diagram, frame, {0, axes.aspect.y, 0}, title + "!!t");
  }

  auto min    = vec3f{axes.xbounds.x, axes.ybounds.x, -1};
  auto max    = vec3f{axes.xbounds.y, axes.ybounds.y, +1};
  auto center = (max + min) / 2, size = (max - min) / 2;
  auto pframe = scaling_frame(vec3f{axes.aspect.x, axes.aspect.y, 1} / size) *
                translation_frame(-center);

  auto ticks = label_data{};
  for (auto idx : range(axes.xlabels.size())) {
    add_label(diagram, frame,
        {transform_point(pframe, {axes.xticks[idx], 0, 0}).x, -axes.aspect.y,
            0},
        axes.xlabels[idx] + "!!b");
  }
  for (auto idx : range(axes.ylabels.size())) {
    add_label(diagram, frame,
        {-axes.aspect.x, transform_point(pframe, {0, axes.yticks[idx], 0}).y,
            0},
        axes.ylabels[idx] + "!!l");
  }

  return frame * pframe;
}

void add_lineplot(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& points, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto clip = invalidb2f;  // TODO: save on context
  if (clip == invalidb2f) {
    auto positions = vector<vec3f>{};
    for (auto point : points) positions.push_back({point.x, point.y, 0});
    return add_polyline(diagram, frame, positions, {}, stroke, thickness);
  } else {
    auto clamp = [](vec2f x, vec2f min, vec2f max) -> vec2f {
      return {yocto::clamp(x.x, min.x, max.x), yocto::clamp(x.y, min.y, max.y)};
    };
    auto positions = vector<vec3f>{};
    for (auto idx : range((int)points.size() - 1)) {
      auto point1 = points[idx + 0], point2 = points[idx + 1];
      auto clipped1 = clamp(point1, clip.min, clip.max),
           clipped2 = clamp(point2, clip.min, clip.max);
      if (point1 != clipped1 && point2 != clipped2) continue;
      positions.push_back({clipped1.x, clipped1.y, 0});
      positions.push_back({clipped2.x, clipped2.y, 0});
    }
    return add_lines(diagram, frame, positions, labels, stroke, thickness);
  }
};
void add_lineplot(diagram_data& diagram, const frame3f& frame,
    const diagram_func<float, float>& function, const vec2f& range,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_lineplot(diagram, frame, sample_function(function, range, 100),
      labels, stroke, thickness);
}

void add_scatterplot(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& points, const vector<string>& labels,
    const vec4f& stroke, float thickness) {
  auto clip = invalidb2f;  // TODO: save on context
  if (clip == invalidb2f) {
    auto positions = vector<vec3f>{};
    for (auto point : points) positions.push_back({point.x, point.y, 0});
    return add_points(diagram, frame, positions, {}, stroke, thickness);
  } else {
    auto clamp = [](vec2f x, vec2f min, vec2f max) -> vec2f {
      return {yocto::clamp(x.x, min.x, max.x), yocto::clamp(x.y, min.y, max.y)};
    };
    auto positions = vector<vec3f>{};
    for (auto point : points) {
      auto clipped = clamp(point, clip.min, clip.max);
      if (point == clipped) positions.push_back({point.x, point.y, 0});
    }
    return add_points(diagram, frame, positions, labels, stroke, thickness);
  }
}

void add_scatterplot(diagram_data& diagram, const frame3f& frame,
    const diagram_func<float, float>& function, const vec2f& range,
    const vector<string>& labels, const vec4f& stroke, float thickness) {
  return add_scatterplot(diagram, frame, sample_function(function, range, 100),
      labels, stroke, thickness);
}

vector<vec2f> sample_function(const diagram_func<float, float>& function,
    const vec2f& range_, int samples) {
  auto data = vector<vec2f>(samples);
  for (auto idx : range(samples)) {
    auto x    = lerp(range_.x, range_.y, (float)idx / (float)(samples - 1));
    data[idx] = {x, function(x)};
  }
  return data;
}

// Color shapes
array2d<vec4f> random_image(const vec2i& steps) {
  auto image = array2d<vec4f>(steps);
  auto rng   = make_rng(17);
  for (auto ij : range(steps)) {
    auto color = rand3f(rng) * 0.6f + 0.4f;
    image[ij]  = {color.x, color.y, color.z, 1};
  }
  return image;
}
array2d<vec4f> gamma_image(const vec2i& steps) {
  auto image = array2d<vec4f>(steps);
  for (auto idx : range((size_t)steps.x)) {
    auto color      = idx / (float)steps.x;
    image[{idx, 0}] = {color, color, color, 1};
  }
  for (auto idx : range((size_t)steps.x)) {
    auto color      = pow(idx / (float)steps.x, 2.2f);
    image[{idx, 1}] = {color, color, color, 1};
  }
  return image;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL API FOR DIAGRAM RENDERING
// -----------------------------------------------------------------------------
namespace yocto {
// crop image vertically
static array2d<vec4f> crop_image(const array2d<vec4f>& image) {
  // find min and max
  auto [width, height] = (vec2s)image.extents();
  auto min = (size_t)0, max = height;
  for (auto j : range(height)) {
    auto empty = true;
    for (auto i : range(width)) {
      if (image[{i, j}] != vec4f{1, 1, 1, 1}) empty = false;
    }
    if (empty) {
      min = j;
    } else {
      break;
    }
  }
  for (auto j : range(height)) {
    auto empty = true;
    for (auto i : range(width)) {
      if (image[{i, height - j - 1}] != vec4f{1, 1, 1, 1}) empty = false;
    }
    if (empty) {
      max = height - j;
    } else {
      break;
    }
  }

  // crop
  auto cropped = array2d<vec4f>({width, max - min});
  for (auto j : range(min, max)) {
    for (auto i : range(width)) {
      cropped[{i, j - min}] = image[{i, j}];
    }
  }

  // done
  return cropped;
}

// Render a diagram to an image
array2d<vec4f> render_diagram(
    const diagram_data& diagram, bool draft, bool crop) {
  // convert diagram to scene
  auto& scene = diagram.scene;

  // rendering params
  auto params       = trace_params{};
  params.resolution = 1440;
  params.samples    = draft ? 16 : 256;
  params.sampler    = trace_sampler_type::diagram;

  // build bvh
  auto bvh = make_trace_bvh(scene, params);

  // init renderer
  auto lights = make_trace_lights(scene, params);

  // state
  auto state = make_trace_state(scene, params);

  // switch between interactive and offline
  // rendering
  for (auto sample : range(0, params.samples, params.batch)) {
    trace_samples(state, scene, bvh, lights, params);
  }

  // image
  auto image = get_image(state);

  // crop
  if (crop) image = crop_image(image);

  // done
  return image;
}

array2d<vec4f> render_diagram(
    const diagram_callback& callback, bool draft, bool crop) {
  auto diagram = diagram_data{};
  auto context = diagram_context{diagram, diagram.scene};
  callback(context);
  return render_diagram(diagram, draft, crop);
}

}  // namespace yocto

#ifdef YOCTO_OPENGL

#include <yocto/yocto_gui.h>

// -----------------------------------------------------------------------------
// HIGH-LEVEL API FOR DIAGRAM RENDERING
// -----------------------------------------------------------------------------
namespace yocto {

// Render a diagram interactively
void view_diagram(const string& name, const diagram_data& diagram) {
  // convert diagram to scene
  auto scene = diagram.scene;

  // rendering params
  auto params       = trace_params{};
  params.resolution = 1440;
  params.sampler    = trace_sampler_type::diagram;

  // render
  show_trace_gui("ydiagram", name, scene, params, 0.5f);
}

void view_diagram(const string& name, const diagram_callback& callback) {
  auto diagram = diagram_data{};
  auto context = diagram_context{diagram, diagram.scene};
  callback(context);
  return view_diagram(name, diagram);
}

#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide catmullclark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_impl_(
    const vector<vec4i>& quads, const vector<T>& vert, bool uncorrected) {
  // get edges
  auto emap     = make_edge_map(quads);
  auto edges    = get_edges(emap);
  auto boundary = get_boundary(emap);
  // number of elements
  auto nverts = (int)vert.size();
  auto nedges = (int)edges.size();
  auto nfaces = (int)quads.size();

  // split elements ------------------------------------
  // create vertices
  auto tvert = vector<T>(nverts + nedges + nfaces);
  for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
  for (auto i = 0; i < nedges; i++) {
    auto e            = edges[i];
    tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
  }
  for (auto i = 0; i < nfaces; i++) {
    auto q = quads[i];
    if (q.z != q.w) {
      tvert[nverts + nedges + i] =
          (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4;
    } else {
      tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.z]) / 3;
    }
  }
  // create quads
  auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
  auto qi     = 0;
  for (auto i = 0; i < nfaces; i++) {
    auto q = quads[i];
    if (q.z != q.w) {
      tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
          nverts + nedges + i, nverts + edge_index(emap, {q.w, q.x})};
      tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
          nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
      tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.w}),
          nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
      tquads[qi++] = {q.w, nverts + edge_index(emap, {q.w, q.x}),
          nverts + nedges + i, nverts + edge_index(emap, {q.z, q.w})};
    } else {
      tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
          nverts + nedges + i, nverts + edge_index(emap, {q.z, q.x})};
      tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
          nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
      tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.x}),
          nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
    }
  }
  tquads.resize(qi);

  // averaging pass ----------------------------------
  auto avert  = vector<T>(tvert.size(), T{});
  auto acount = vector<int>(tvert.size(), 0);
  for (auto& q : tquads) {
    auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4;
    for (auto vid : {q.x, q.y, q.z, q.w}) {
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];

  // correction pass ----------------------------------
  // p = p + (avg_p - p) * (4/avg_count)
  if (!uncorrected) {
    for (auto i = 0; i < tvert.size(); i++) {
      avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4 / (float)acount[i]);
    }
  }
  tvert = avert;

  return {tquads, tvert};
}

// Bspline subdivision
pair<vector<vec4i>, vector<float>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<float>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec2f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec3f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec4f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}

// Subdivide bspline.
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_bspline_impl(
    const vector<vec2i>& lines, const vector<T>& vert) {
  // split elements ------------------------------------
  auto [tlines, tvert] = subdivide_lines(lines, vert);

  // boundary ------------------------------------------
  auto valence = vector<bool>(tvert.size(), 1);
  // TODO: fixme
  auto tcreases = vector<int>{};

  // averaging pass ----------------------------------
  auto avert  = vector<T>(tvert.size(), T{});
  auto acount = vector<int>(tvert.size(), 0);
  for (auto& p : tcreases) {
    auto c = tvert[p];
    for (auto vid : {p}) {
      if (valence[vid] != 0) continue;
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto& l : tlines) {
    auto c = (tvert[l.x] + tvert[l.y]) / 2;
    for (auto vid : {l.x, l.y}) {
      if (valence[vid] != 1) continue;
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];
  tvert = avert;

  return {tlines, tvert};
}

pair<vector<vec2i>, vector<float>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<float>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec2f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec2f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec3f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec3f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec4f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec4f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE DIAGRAMS
// -----------------------------------------------------------------------------
namespace yocto {

using namespace dcolors;
using namespace dthickness;
using namespace dconstants;

diagram_data placeholder_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"missing"})) {
    auto group = add_group(diagram);
    add_quad(diagram, group, {-2.5, 0.0, 0.0}, 1, {"missing!!c"}, fill1);
    add_quad(diagram, group, {0, 0, 0}, 1, {"missing!!c"}, fill2);
    add_quad(diagram, group, {2.5, 0.0, 0.0}, 1, {"missing!!c"}, fill3);
  }

  return diagram;
}

diagram_data frame_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"composition"})) {
    auto group1 = add_group(diagram, {0.0, 0.0, 0.0});
    add_point(diagram, group1, {0.0, 0.77, 0.0},
        "$\\mathbf{p}_1 \\equiv \\mathbf{p}_2$!!t");
    add_label(diagram, group1, {0.0, 0.77, 0.0},
        "$p_{1,x} \\neq p_{2,x}$ \\\\"
        "$p_{1,y} \\neq p_{2,y}$ \\\\"
        "$p_{1,z} \\neq p_{2,z}$!!b");
    auto group2 = add_group(diagram, {-1.5, 0.0, 0.0});
    add_label(diagram, group2, {0.0, -1.0, 0.0},
        "$\\mathbf{F}_1 = (\\mathbf{x}_1,\\mathbf{y}_1,\\mathbf{z}_1,\\mathbf{o}_1)$");
    add_axes(diagram, group2, drotation({30.0, 205, 0.0}), 1.5,
        {"$\\mathbf{o}_1$!!b", "$\\mathbf{x}_1$!!t", "$\\mathbf{y}_1$!!r",
            "$\\mathbf{z}_1$!!l"},
        stroke1);
    auto group3 = add_group(diagram, {+1.5, 0.0, 0.0});
    add_label(diagram, group3, {0.0, -1.0, 0.0},
        "$\\mathbf{F}_2 = (\\mathbf{x}_2,\\mathbf{y}_2,\\mathbf{z}_2,\\mathbf{o}_2)$");
    add_axes(diagram, group3, drotation({30.0, 105, 25}), 1.5,
        {"$\\mathbf{o}_2$!!b", "$\\mathbf{x}_2$!!r", "$\\mathbf{y}_2$!!r",
            "$\\mathbf{z}_2$!!t"},
        stroke3);
  }
  if (contains(tag, {"coordinates"})) {
    auto group = add_group(diagram, {0, 0, 0}, drotation({30.0, 300.0, 0.0}));
    add_axes(diagram, group, {0, 0, 0}, 1.5,
        {"$\\mathbf{o}$!!b", "$\\mathbf{x}$!!r", "$\\mathbf{y}$!!r",
            "$\\mathbf{z}$!!l"},
        black);
    add_point(diagram, group, {1.0, 1.0, 1.0}, "$\\mathbf{p}$!!t");
    add_label(diagram, group, {1.0, 0.0, 0.0}, "$p_x$!!tr");
    add_label(diagram, group, {0.0, 1.0, 0.0}, "$p_y$!!tr");
    add_label(diagram, group, {0.0, 0.0, 1.0}, "$p_z$!!tl");
    add_cube(diagram, group, {0.5, 0.5, 0.5}, 0.5, {}, transparent, stroke1);
  }
  if (contains(tag, {"definition"})) {
    auto group1 = add_group(diagram, {0.0, 0.0, 0.0});
    add_axes(diagram, group1, drotation({30.0, 300.0, 0.0}), 1,
        {"$\\mathbf{o}$!!b", "$\\mathbf{x}$!!r", "$\\mathbf{y}$!!r",
            "$\\mathbf{z}$!!l"},
        black);
    auto group2 = add_group(diagram, {-1.75, 0.0, 0.0});
    add_label(diagram, group2, {0.0, -1.0, 0.0},
        "$\\mathbf{F}_1 = (\\mathbf{x}_1,\\mathbf{y}_1,\\mathbf{z}_1,\\mathbf{o}_1)$");
    add_axes(diagram, group2, drotation({30.0, 205, 0.0}), 1,
        {"$\\mathbf{o}_1$!!b", "$\\mathbf{x}_1$!!t", "$\\mathbf{y}_1$!!r",
            "$\\mathbf{z}_1$!!l"},
        stroke1);
    auto group3 = add_group(diagram, {+1.75, 0.0, 0.0});
    add_label(diagram, group3, {0.0, -1.0, 0.0},
        "$\\mathbf{F}_2 = (\\mathbf{x}_2,\\mathbf{y}_2,\\mathbf{z}_2,\\mathbf{o}_2)$");
    add_axes(diagram, group3, drotation({30.0, 105, 25}), 1,
        {"$\\mathbf{o}_2$!!b", "$\\mathbf{x}_2$!!r", "$\\mathbf{y}_2$!!r",
            "$\\mathbf{z}_2$!!t"},
        stroke2);
  }
  return diagram;
}

diagram_data objtransform_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"rotationx"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1, {}, black);
      add_axes(diagram, group1, drotation({30.0, 0.0, 0.0}), 1, {}, stroke3);
      add_arc(diagram, group1 * drotation({0.0, 90.0, 90.0}), {0, 0, 0}, 1,
          radians(45), {"$\\theta$!!ltc"}, stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (s \\mathbf{x}, s \\mathbf{y}, s \\mathbf{z}, \\mathbf{o} ) $");
    }
    {
      auto group3 = add_group(
          diagram, {0, 0, 0}, dtransform({2.5, 0.0, 0.0}, {15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3 * drotation({30.0, 0.0, 0.0}), {0, 0, 0}, 0.5);
      add_axes(diagram, group3 * drotation({30.0, 0.0, 0.0}), {0, 0, 0}, 1, {},
          stroke3);
    }
  }
  if (contains(tag, {"rotationy"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1, drotation({0.0, 60.0, 0.0}), 1, {}, stroke3);
      add_arc(diagram, group1 * drotation({270.0, 0.0, 0.0}), {0, 0, 0}, 1,
          radians(45), {"$\\theta$!!rrrc"}, stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (s \\mathbf{x}, s \\mathbf{y}, s \\mathbf{z}, \\mathbf{o} ) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
    }
    {
      auto group3 = add_group(
          diagram, {2.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3 * drotation({0.0, 60.0, 0.0}), {0, 0, 0}, 0.5);
      add_axes(diagram, group3 * drotation({0.0, 60.0, 0.0}), {0, 0, 0}, 1, {},
          stroke3);
    }
  }
  if (contains(tag, {"rotationz"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1 * drotation({0.0, 0.0, 60.0}), {0, 0, 0}, 1, {},
          stroke3);
      add_arc(diagram, group1, {0, 0, 0}, 1, radians(45), {"$\\theta$!!rttc"},
          stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (s \\mathbf{x}, s \\mathbf{y}, s \\mathbf{z}, \\mathbf{o} ) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
    }
    {
      auto group3 = add_group(
          diagram, {2.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3 * drotation({0.0, 0.0, 60.0}), {0, 0, 0}, 0.5);
      add_axes(diagram, group3 * drotation({0.0, 0.0, 60.0}), {0, 0, 0}, 1, {},
          stroke3);
    }
  }
  if (contains(tag, {"scalingnu"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1 * dscaling({0.9, 0.6, 0.8}), {0, 0, 0}, 1, {},
          stroke3);
      add_points(diagram, group1,
          {{0.9, 0.0, 0.0}, {0.0, 0.6, 0.0}, {0.0, 0.0, 0.8}},
          {"$s_x$!!b", "$s_y$!!l", "$s_z$!!b"}, stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (s_x \\mathbf{x}, s_y \\mathbf{y}, s_z \\mathbf{z}, \\mathbf{o} ) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
    }
    {
      auto group3 = add_group(
          diagram, {2.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3 * dscaling({0.9, 0.6, 0.8}), {0, 0, 0}, 0.5);
      add_axes(diagram, group3 * dscaling({0.9, 0.6, 0.8}), {0, 0, 0}, 1, {},
          stroke3);
    }
  }
  if (contains(tag, {"scalingu"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1 * dscaling({0.8, 0.8, 0.8}), {0, 0, 0}, 1, {},
          stroke3);
      add_points(diagram, group1,
          {{0.8, 0.0, 0.0}, {0.0, 0.8, 0.0}, {0.0, 0.0, 0.8}},
          {"$s$!!b", "$s$!!l", "$s$!!b"}, stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (s \\mathbf{x}, s \\mathbf{y}, s \\mathbf{z}, \\mathbf{o} ) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
    }
    {
      auto group3 = add_group(
          diagram, {2.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3 * dscaling({0.8, 0.8, 0.8}), {0, 0, 0}, 0.5);
      add_axes(diagram, group3 * dscaling({0.8, 0.8, 0.8}), {0, 0, 0}, 1, {},
          stroke3);
    }
  }
  if (contains(tag, {"translation"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1, {0.65, 0.0, 0.3}, 1, {}, stroke3);
      add_arrow(diagram, group1, {0.0, 0.0, 0.0}, {0.65, 0.0, 0.3},
          {"$\\mathbf{t}$!!le"}, stroke2);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$ \\mathbf{F} = (\\mathbf{x}, \\mathbf{y}, \\mathbf{z}, \\mathbf{o} + \\mathbf{t} ) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_cube(diagram, group2, {0, 0, 0}, 0.5);
    }
    {
      auto group3 = add_group(
          diagram, {2.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_cube(diagram, group3, {0.65, 0.0, 0.3}, 0.5);
      add_axes(diagram, group3, {0.65, 0.0, 0.3}, 1, {}, stroke3);
    }
  }
  return diagram;
}

diagram_data primitive_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"types"})) {
    {
      auto group1 = add_group(diagram, {-2.5, 0.0, 0.0});
      add_point(diagram, group1, {0.0, 0.0, 0.0}, "$\\mathbf{p}_1$!!r");
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "point");
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_arrow(diagram, group2, {-0.35, -0.7, 0.0}, {0.35, 0.7, 0.0},
          {"$\\mathbf{v}$!!rc"}, black);
      add_points(diagram, group2, {{-0.35, -0.7, -0.03}, {0.35, 0.7, -0.03}},
          {"$\\mathbf{p}_1$!!r", "$\\mathbf{p}_2$!!r"}, stroke1);
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "vector");
    }
    {
      auto group3 = add_group(diagram, {2.5, 0.0, 0.0});
      add_arrow(diagram, group3, {-0.35, -0.7, 0.0}, {0.15, 0.3, 0.0},
          {"$\\mathbf{d}$!!rc"}, black);
      add_arrow(diagram, group3, {-0.35, -0.7, -0.01}, {0.35, 0.7, -0.01},
          {"$\\mathbf{v}$!!rrrrttttttttc"}, stroke1);
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "direction");
    }
  }
  return diagram;
}

diagram_data transform_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"rotation"})) {
    {
      auto group1 = add_group(
          diagram, {-2.0, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_label(diagram, group1, {0.0, -1.25, 0.0},
          "$ \\mathbf{F} = (c \\mathbf{x} + s \\mathbf{y}, s \\mathbf{y} + c \\mathbf{z}, 0) $!!l");
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_axes(diagram, group1 * drotation({45.0, 0.0, 0.0}), {0, 0, 0}, 1, {},
          stroke1);
      add_arc(diagram, group1 * drotation({315.0, 270.0, 0.0}), {0, 0, 0}, 0.5,
          radians(45), {}, stroke3);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -1.25, 0.0},
          "$ \\mathbf{F} = (c \\mathbf{x} + s \\mathbf{y}, s \\mathbf{y} + c \\mathbf{z}, 0) $");
      add_axes(diagram, group2, {0, 0, 0}, 1);
      add_axes(diagram, group2 * drotation({0.0, 45.0, 0.0}), {0, 0, 0}, 1, {},
          stroke1);
      add_arc(diagram, group2 * drotation({90.0, 0.0, 45.0}), {0, 0, 0}, 0.5,
          radians(45), {}, stroke3);
    }
    {
      auto group3 = add_group(
          diagram, {2.0, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group3, {0.0, -1.25, 0.0},
          "$ \\mathbf{F} = (c \\mathbf{x} + s \\mathbf{y}, s \\mathbf{y} + c \\mathbf{z}, 0) $");
      add_axes(diagram, group3, {0, 0, 0}, 1);
      add_axes(diagram, group3 * drotation({0.0, 0.0, 45.0}), {0, 0, 0}, 1, {},
          stroke1);
      add_arc(diagram, group3, {0, 0, 0}, 0.5, radians(45), {}, stroke3);
    }
  }
  if (contains(tag, {"scaling"})) {
    {
      auto group1 = add_group(
          diagram, {-1.5, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_label(diagram, group1, {0.0, -0.77, 0.0},
          "$ \\mathbf{F} = (s \\mathbf{x}, s \\mathbf{y}, s \\mathbf{z} ) $");
      add_axes(diagram, group1, {0, 0, 0}, 1);
      add_lines(diagram, group1,
          {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.5, 0.0},
              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.5}});
      add_points(diagram, group1,
          {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}}, {}, stroke3);
    }
    {
      auto group2 = add_group(
          diagram, {1.5, 0.0, 0.0}, drotation({15.0, 325.0, 0.0}));
      add_label(diagram, group2, {0.0, -0.77, 0.0},
          "$ \\mathbf{F} = (s_x \\mathbf{x}, s_y \\mathbf{y}, s_z \\mathbf{z} ) $");
      add_lines(diagram, group2,
          {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.75, 0.0},
              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.5}},
          {}, stroke1);
      add_points(diagram, group2,
          {{0.25, 0.0, 0.0}, {0.0, 0.75, 0.0}, {0.0, 0.0, 0.5}}, {}, stroke3);
      add_axes(diagram, group2, {0, 0, 0}, 1);
    }
  }
  if (contains(tag, {"translation"})) {
    {
      auto group = add_group(
          diagram, {-1.5, 0.0, 0.0}, drotation({30.0, 315.0, 0.0}));
      add_label(diagram, group, {0.0, -0.77, 0.0},
          "$ \\mathbf{F} = (\\mathbf{x}, \\mathbf{y}, \\mathbf{z}, \\mathbf{o} + \\mathbf{t} ) $");
      add_axes(diagram, group, {0, 0, 0}, 1);
      add_axes(diagram, group * dtranslation({0.385, 0.0, -0.75}), {0, 0, 0}, 1,
          {}, stroke1);
      add_arrow(diagram, group, {-1.5, 0.0, 0.0}, {-0.8, 0.1, 0.0},
          {"$\\mathbf{t}$!!tc"}, stroke3);
    }
  }
  return diagram;
}

diagram_data vector_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"operations"})) {
    {
      auto group1 = add_group(diagram, {-3.75, 0.0, 0.0});
      add_vector(
          diagram, group1, {1.5, 0.5, 0.0}, {"$\\mathbf{u}+\\mathbf{v}$!!bre"});
      add_vector(
          diagram, group1, {0.5, 1.0, 0.0}, {"$\\mathbf{v}$!!re"}, stroke1);
      add_vector(
          diagram, group1, {1.0, -0.5, 0.0}, {"$\\mathbf{u}$!!ble"}, stroke2);
    }
    {
      auto group2 = add_group(diagram, {-1.0, 0.0, 0.0});
      add_vector(diagram, group2, {0.25, 0.5, 0.0}, {"$s\\mathbf{v}$!!re"});
      add_vector(
          diagram, group2, {0.5, 1.0, 0.0}, {"$\\mathbf{v}$!!trre"}, stroke1);
    }
    {
      auto group3 = add_group(diagram, {1.0, 0.0, 0.0});
      add_vector(diagram, group3, {-0.5, -1.0, 0.0}, {"$-\\mathbf{v}$!!re"});
      add_vector(
          diagram, group3, {0.5, 1.0, 0.0}, {"$\\mathbf{v}$!!re"}, stroke1);
    }
    {
      auto group4 = add_group(diagram, {3.0, 0.0, 0.0});
      add_vector(diagram, group4, {-0.5, 0.5, 0.0},
          {"$\\mathbf{v}-\\mathbf{u}$!!ble"});
      add_vector(
          diagram, group4, {0.5, 1.0, 0.0}, {"$\\mathbf{v}$!!trre"}, stroke1);
      add_vector(
          diagram, group4, {1.0, 0.5, 0.0}, {"$\\mathbf{u}$!!be"}, stroke2);
    }
  }
  if (contains(tag, {"products"})) {
    {
      auto group1  = add_group(diagram, {-1.5, 0.0, 0.0});
      auto rgroup1 = group1 * drotation({50.0, 15.0, 315.0});
      add_label(diagram, group1, {0.0, -1.0, 0.0},
          "$| \\mathbf{u} \\times \\mathbf{v} | = |\\mathbf{u}| |\\mathbf{v}| \\sin\\theta$!!");
      add_vector(
          diagram, rgroup1, {0.0, 1.25, 0.0}, {"$\\mathbf{v}$!!tle"}, stroke1);
      add_vector(
          diagram, rgroup1, {1.25, 0.0, 0.0}, {"$\\mathbf{u}$!!ble"}, stroke2);
      add_arc(diagram, rgroup1, {0, 0, 0}, 0.5, radians(90), {}, stroke3);
      add_label(diagram, rgroup1, {0.5, 0.5, 0.0}, "$\\theta$");
    }
    {
      auto group2  = add_group(diagram, {1.0, 0.0, 0.0});
      auto rgroup2 = group2 * drotation({50.0, 15.0, 315.0});
      add_label(diagram, group2, {0.0, -1.0, 0.0},
          "$\\mathbf{u} \\cdot \\mathbf{v} = |\\mathbf{u}| |\\mathbf{v}| \\cos\\theta$!!");
      add_vector(diagram, rgroup2, {0.0, 0.0, -1.5},
          {"$\\mathbf{u} \\times \\mathbf{v}$!!le"});
      add_vector(
          diagram, rgroup2, {0.0, 1.25, 0.0}, {"$\\mathbf{v}$!!tle"}, stroke1);
      add_vector(
          diagram, rgroup2, {1.25, 0.0, 0.0}, {"$\\mathbf{u}$!!ble"}, stroke2);
      add_arc(diagram, rgroup2, {0, 0, 0}, 0.5, radians(90), {}, stroke3);
      add_label(diagram, rgroup2, {0.5, 0.5, 0.0}, "$\\theta$");
    }
  }
  return diagram;
}

diagram_data compositing_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto bottom = vector<vec3f>{
      {-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
  auto top = vector<vec3f>{
      {-1.0, 0.2, 0.0}, {1.0, -0.2, 0.0}, {1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
  auto top_reminder = vector<vec3f>{
      {-1.0, 0.2, 0.0}, {1.0, -0.2, 0.0}, {1.0, -1.0, 0.0}, {-1.0, -1.0, 0.0}};
  auto mid = vector<vec3f>{
      {-1.0, -1.0, 0.0}, {-0.4, -1.0, 0.0}, {0.4, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
  auto mid_reminder = vector<vec3f>{
      {1.0, -1.0, 0.0}, {-0.4, -1.0, 0.0}, {0.4, 1.0, 0.0}, {1.0, 1.0, 0.0}};
  auto topmid_reminder = vector<vec3f>{
      {1.0, -1.0, 0.0}, {-0.4, -1.0, 0.0}, {0.4, 1.0, 0.0}, {1.0, 1.0, 0.0}};

  if (contains(tag, {"alpha1"})) {
    {
      auto group1 = add_group(diagram, {-3.3, 0.0, 0.0});
      add_quadv(diagram, group1, top, {"$\\alpha_A$!!e"}, fill1);
      add_quadv(diagram, group1, top_reminder, {"$1-\\alpha_A$!!e"}, fill3);
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "$A$!!t");
    }
    {
      auto group2 = add_group(diagram, {-1.1, 0.0, 0.0});
      add_quadv(diagram, group2, bottom, {"$\\alpha_B=1$!!e"}, fill3);
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "$B$!!t");
    }
    {
      auto group3 = add_group(diagram, {1.1, 0.0, 0.0});
      add_quadv(diagram, group3, top, {"$\\alpha_A$!!e"}, fill1);
      add_quadv(diagram, group3, top_reminder, {"$1-\\alpha_A$!!e"}, fill3);
      add_label(
          diagram, group3, {0.0, 1.0, 0.0}, "$A \\ \\text{over} \\ B$!!t");
    }
    {
      auto group4 = add_group(
          diagram, {3.3, 0.0, 0.0}, drotation({21.0, 332.0, 0.0}));
      add_quadv(diagram, group4, translate({0.0, 0.0, 0.4}, top), {}, fill1);
      add_quadv(diagram, group4, bottom, {}, fill3);
    }
  }
  if (contains(tag, {"alpha2"})) {
    {
      auto group1 = add_group(diagram, {-3.3, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "$A$!!t");
      add_quadv(diagram, group1, top, {"$\\alpha_A$!!e"}, fill1);
      add_quadv(diagram, group1, top_reminder, {"$1-\\alpha_A$!!e"}, fill3);
    }
    {
      auto group2 = add_group(diagram, {-1.1, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "$B$!!t");
      add_quadv(diagram, group2, mid, {"$\\alpha_B$!!e"}, fill2);
      add_quadv(diagram, group2, mid_reminder, {"$1-\\alpha_B$!!e"}, fill3);
    }
    {
      auto group3 = add_group(diagram, {1.1, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0},
          "$A \\ \\text{over} \\ B \\ \\text{over} \\ C$!!t");
      add_quadv(diagram, group3, top, {"$\\alpha_A$!!e"}, fill1);
      add_quadv(diagram, group3,
          {{0.0, 0.0, 0.0}, {1.0, -0.2, 0.0}, {1.0, -1.0, 0.0},
              {-0.4, -1.0, 0.0}},
          {"$(1-\\alpha_A) \\cdot$ \\\\ $(1-\\alpha_B)$!!e"}, fill3);
      add_quadv(diagram, group3,
          {{-1.0, -1.0, 0.0}, {-0.4, -1.0, 0.0}, {0.0, 0.0, 0.0},
              {-1.0, 0.2, 0.0}},
          {"$(1-\\alpha_A) \\cdot$ \\\\ $\\alpha_B$!!e"}, fill2);
    }
    {
      auto group4 = add_group(
          diagram, {3.3, 0.0, 0.0}, drotation({21.0, 332.0, 0.0}));
      add_quadv(diagram, group4, translate({0.0, 0.0, 0.4}, top), {}, fill1);
      add_quadv(diagram, group4, bottom, {}, fill3);
      add_quadv(diagram, group4, translate({0.0, 0.0, 0.2}, mid), {}, fill2);
    }
  }
  if (contains(tag, {"alpha3"})) {
    {
      auto group1 = add_group(diagram, {-3.3, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "coverage!!t");
      add_quadv(diagram, group1, top, {}, fill1);
      add_quadv(diagram, group1, top_reminder, {}, fill3);
    }
    {
      auto group2 = add_group(
          diagram, {-1.1, 0.0, 0.0}, drotation({21, 335, 0.0}));
      add_quadv(diagram, group2, translate({0.0, 0.0, 0.4}, top), {}, fill1);
      add_quadv(diagram, group2, bottom, {}, fill3);
    }
    {
      auto group3 = add_group(diagram, {1.1, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "opacity!!t");
      add_quadv(
          diagram, group3, translate({0.0, 0.0, 0.02}, bottom), {}, tfill1);
      add_quadv(diagram, group3, bottom, {}, fill3);
    }
    {
      auto group4 = add_group(
          diagram, {3.3, 0.0, 0.0}, drotation({21.0, 332.0, 0.0}));
      add_quadv(
          diagram, group4, translate({0.0, 0.0, 0.4}, bottom), {}, tfill1);
      add_quadv(diagram, group4, bottom, {}, fill3);
    }
  }
  if (contains(tag, {"alpha4"})) {
    {
      auto group1 = add_group(diagram, {-3.3, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "valid over!!t");
      add_quadv(diagram, group1, top, {}, fill1);
      add_quadv(diagram, group1,
          {{0.0, 0.0, 0.0}, {1.0, -0.2, 0.0}, {1.0, -1.0, 0.0},
              {-0.4, -1.0, 0.0}},
          {}, fill3);
      add_quadv(diagram, group1,
          {{-1.0, -1.0, 0.0}, {-0.4, -1.0, 0.0}, {0.0, 0.0, 0.0},
              {-1.0, 0.2, 0.0}},
          {}, fill2);
    }
    {
      auto group2 = add_group(
          diagram, {-1.1, 0.0, 0.0}, drotation({21, 335, 0.0}));
      add_quadv(diagram, group2, translate({0.0, 0.0, 0.2}, mid), {}, fill2);
      add_quadv(diagram, group2, translate({0.0, 0.0, 0.4}, top), {}, fill1);
      add_quadv(diagram, group2, bottom, {}, fill3);
    }
    {
      auto group3 = add_group(diagram, {1.1, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "invalid over!!t");
      add_quadv(diagram, group3, top, {}, fill1);
      add_quadv(diagram, group3, top, {}, fill1);
      add_quadv(diagram, group3,
          {{-1.0, 0.2, 0.0}, {1.0, -0.2, 0.0}, {1.0, -0.6, 0.0},
              {-1.0, -0.2, 0.0}},
          {}, fill2);
      add_quadv(diagram, group3,
          {{-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, -0.6, 0.0},
              {-1.0, -0.2, 0.0}},
          {}, fill3);
    }
    {
      auto group4 = add_group(
          diagram, {3.3, 0.0, 0.0}, drotation({21.0, 332.0, 0.0}));
      add_quadv(diagram, group4, translate({0.0, 0.0, 0.4}, top), {}, fill1);
      add_quadv(diagram, group4, bottom, {}, fill3);
      add_quadv(diagram, group4,
          translate({0.0, 0.0, 0.2}, {{-1.0, -0.2, 0.0}, {1.0, -0.6, 0.0},
                                         {1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}}),
          {}, fill2);
    }
  }
  return diagram;
}

diagram_data image_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto image = load_image(
      "tests/_diagramdata/images/unsplash/eberhard-grossgasteiger-unsplash.jpg");

  if (contains(tag, {"approximation"})) {
    {
      // TODO: image borders
      auto group1 = add_group(diagram, {-1.93, 0.0, 0.0});
      add_image(diagram, group1, {0, 0, 0}, 1.2, image);
      add_label(diagram, group1, {0.0, 1.2, 0.0}, "image!!t");
    }
    {
      // TODO: image borders
      auto group2 = add_group(diagram, {1.93, 0.0, 0.0});
      add_image(diagram, group2, {0, 0, 0}, 1.2, resize_image(image, {32, 0}));
      add_label(diagram, group2, {0.0, 1.2, 0.0}, "sampled image!!t");
    }
  }
  if (contains(tag, {"grid"})) {
    {
      auto group1 = add_group(diagram, {-2.1, 0.0, 0.0});
      add_imagegrid(diagram, group1, {0, 0, 0}, 1, random_image({4, 2}),
          {"$(0,0)$", "$(1,0)$", "$(2,0)$", "$(3,0)$", "$(0,1)$", "$(1,1)$",
              "$(2,1)$", "$(3,1)$"});
      add_arrows(diagram, group1,
          {{-2.15, 1.15, 0.0}, {2.15, 1.15, 0.0}, {-2.15, 1.15, 0.0},
              {-2.15, -1.15, 0.0}},
          {"", "$i$!!t", "", "$j$!!l"});
    }
    {
      auto group2 = add_group(diagram, {2.3, 0.0, 0.0});
      add_imagegrid(diagram, group2, {0, 0, 0}, 1, random_image({4, 2}),
          {"$\\mathbf{c}{0}$", "$\\mathbf{c}{1}$", "$\\mathbf{c}{2}$",
              "$\\mathbf{c}{3}$", "$\\mathbf{c}{4}$", "$\\mathbf{c}{5}$",
              "$\\mathbf{c}{6}$", "$\\mathbf{c}{7}$"});
    }
  }
  if (contains(tag, {"gamma"})) {
    auto group = add_group(diagram);
    add_imagegrid(diagram, group, {0, 0, 0}, 0.8, gamma_image({10, 2}));
    add_label(diagram, group, {0.0, 0.8, 0.0}, "gamma-corrected!!t");
    add_label(diagram, group, {0.0, -0.8, 0.0}, "linear!!b");
  }
  return diagram;
}

// falsecolor
static auto falsecolor_image(const array2d<vec4f>& image) {
  auto minmax = vec2f{flt_max, flt_min};
  for (auto& pixel : image) {
    auto l = std::log10(luminance(xyz(pixel)));
    minmax = {min(minmax.x, l), max(minmax.y, l)};
  }
  auto falsecolor = image;
  for (auto& pixel : falsecolor) {
    auto l = std::log10(luminance(xyz(pixel)));
    auto c = colormap(
        (l - minmax.x) / (minmax.y - minmax.x), colormap_type::viridis);
    pixel = {c.x, c.y, c.z, pixel.w};
  }
  return falsecolor;
}

diagram_data tonemapping_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto image = load_image(
      "tests/_diagramdata/images/hdrihaven/aerodynamicsworkshop_hdrihaven.hdr");

  if (contains(tag, {"exposure"})) {
    auto group = add_group(diagram, {0.0, 0.0, 0.0});
    add_imagemosaic(diagram, group,
        {tonemap_image(image, -2), tonemap_image(image, -1),
            tonemap_image(image, +1), tonemap_image(image, +2)},
        {"$-2$", "$-1$", "$+1$", "$+2$"});
  }
  if (contains(tag, {"falsecolor"})) {
    auto group = add_group(diagram);
    add_imagefull(diagram, group, falsecolor_image(image), "$log(L)$");
  }
  if (contains(tag, {"filmic"})) {
    auto group = add_group(diagram);
    add_imagefull(diagram, group, tonemap_image(image, 0, true), "filmic");
  }
  if (contains(tag, {"gamma"})) {
    auto group = add_group(diagram);
    add_imagefull(diagram, group, tonemap_image(image, 0), "gamma");
  }
  if (contains(tag, {"linear"})) {
    // hack for now
    auto group = add_group(diagram);
    add_imagefull(diagram, group, image, "linear");
  }
  if (contains(tag, {"srgb"})) {
    auto group = add_group(diagram);
    add_imagefull(diagram, group, tonemap_image(image, 0), "sRGB");
  }
  return diagram;
}

diagram_data tonemapplot_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"filmic"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"filmic tonemapping"},
             .xbounds = {0.0, 2.0},
             .xlabels = {"0", "2"},
             .xticks  = {0.0, 2.0},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot,
        [](float x) {
          return rgb_to_srgb(tonemap_filmic(vec3f{x, x, x}).x);
        },
        {0.0, 2.0}, {}, stroke1);
  }
  if (contains(tag, {"gamma"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"gamma correction"},
             .xbounds = {0.0, 2.0},
             .xlabels = {"0", "2"},
             .xticks  = {0.0, 2.0},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot,
        [](float x) { return clamp(pow(x, 1 / 2.2f), 0.0f, 1.0f); }, {0.0, 2.0},
        {}, stroke1);
  }
  if (contains(tag, {"srgb"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"sRGB encoding"},
             .xbounds = {0.0, 2.0},
             .xlabels = {"0", "2"},
             .xticks  = {0.0, 2.0},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot,
        [](float x) { return clamp(rgb_to_srgb(x), 0.0f, 1.0f); }, {0.0, 2.0},
        {}, stroke1);
  }
  if (contains(tag, {"srgbvsgamma"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"sRGB encoding vs gamma correction"},
             .xbounds = {0.0, 2.0},
             .xlabels = {"0", "2"},
             .xticks  = {0.0, 2.0},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot,
        [](float x) { return clamp(rgb_to_srgb(x), 0.0f, 1.0f); }, {0.0, 2.0},
        {}, stroke1);
    add_lineplot(
        diagram, plot,
        [](float x) { return clamp(pow(x, 1 / 2.2f), 0.0f, 1.0f); }, {0.0, 2.0},
        {}, stroke1);
  }
  return diagram;
}

diagram_data camera_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"frame"})) {
    auto group = add_group(diagram, {0, 0, 0}, drotation({20.0, 300.0, 0.0}));
    add_axes(diagram, group, {0, 0, 0}, 1,
        {"$\\mathbf{o}$!!bl", "$\\mathbf{x}$!!b", "$\\mathbf{y}$!!l",
            "$\\mathbf{z}$!!t"});
    add_disk(diagram, group, {0, 0, 0}, 0.2, {}, fill3);
    add_label(diagram, group, {0.0, 0.0, 0.0}, "lens!!tll");
    add_rect(diagram, group, {3, 2}, {0.0, 0.0, 1.5}, 0.4,
        {"", "", "", "film!!l"}, tfill1);
    add_rect(diagram, group, {3, 2}, {0.0, 0.0, -3.0}, 0.8,
        {"", "focus plane!!hbhl"}, tfill2);
    add_lines(diagram, group,
        {{-1.5, -1.0, -3.75}, {0.0, 0.0, 0.0}, {1.5, -1.0, -3.75},
            {0.0, 0.0, 0.0}, {1.5, 1.0, -3.75}, {0.0, 0.0, 0.0},
            {-1.5, 1.0, -3.75}, {0.0, 0.0, 0.0}},
        {}, gray);
    add_label(diagram, group, {-0.6, 0.4, -1.85}, "view pyramid!!ttl");
    add_line(diagram, group, {0.0, 0.0, 1.1}, {0.0, 0.0, 1.5}, {}, gray);
  }
  return diagram;
}

[[maybe_unused]] static shape_data datgrid(const vec2i& steps) {
  auto shape = shape_data{};
  for (auto i = -4; i <= 8; i++) {
    auto u = i / (float)8;
    auto a = vec3f{-1 + u * 2 - 0.075f, -1.1f, 0},
         b = vec3f{u * 2 + 0.075f, 1.1f, 0};
    if (a.x < -1.1f) {
      auto t = (-1.1f - a.x) / (b.x - a.x);
      a      = a + t * (b - a);
    }
    if (b.x > +1.1f) {
      auto t = (+1.1f - a.x) / (b.x - a.x);
      b      = a + t * (b - a);
    }
    shape.lines.push_back(
        {(int)shape.positions.size(), (int)shape.positions.size() + 1});
    shape.positions.push_back(a);
    shape.positions.push_back(b);
  }
  for (auto i = 0; i <= 8; i++) {
    auto u = i / (float)8;
    shape.lines.push_back(
        {(int)shape.positions.size(), (int)shape.positions.size() + 1});
    shape.positions.push_back({-1.1f, -1 + u * 2, 0});
    shape.positions.push_back({1.1f, -1 + u * 2, 0});
  }
  return shape;
}

diagram_data barycentric_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"triangle"})) {
    auto vertices = dconstants::triangle_positions;
    auto point    = vec3f{0.0, -0.2, 0.0};
    {
      auto group1 = add_group(diagram, {-3.0, 0.0, 0.0});
      add_trianglev(diagram, group1, vertices,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_{2}$!!r",
              "$\\mathbf{p}_{3}$!!r"});
      add_point(diagram, group1, point, "$\\mathbf{p}(u,v)$!!t");
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "triangle point!!tt");
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_trianglev(diagram, group2, vertices, {}, fill1);
      add_point(diagram, group2, vertices[0], "$\\mathbf{p}_1$!!l");
      add_arrow(diagram, group2, vertices[0], vertices[1],
          {"", "$\\mathbf{e}_{2,1}$!!r"}, black);
      add_arrow(diagram, group2, vertices[0], vertices[2],
          {"", "$\\mathbf{e}_{3,1}$!!r"}, black);
      add_point(diagram, group2, point, "$\\mathbf{p}(u,v)$!!t");
      add_line(diagram, group2, point, {-0.6, -0.2, 0.0}, {"$u$!!bc"}, stroke3);
      add_line(diagram, group2, point, {-0.4, -1.0, 0.0}, {"$v$!!rc"}, stroke3);
      add_label(
          diagram, group2, {0.0, 1.0, 0.0}, "barycentric coordinates!!tt");
    }
    {
      auto group3 = add_group(diagram, {3.0, 0.0, 0.0});
      add_trianglev(diagram, group3, vertices, {}, fill1);
      add_point(diagram, group3, vertices[0], "$\\mathbf{p}_1$!!l");
      add_arrow(diagram, group3, vertices[0], vertices[1],
          {"", "$\\mathbf{e}_{2,1}$!!r"}, black);
      add_arrow(diagram, group3, vertices[0], vertices[2],
          {"", "$\\mathbf{e}_{3,1}$!!r"}, black);
      add_point(diagram, group3, point, "$\\mathbf{p}(u,v)$!!t");
      add_atgrid(diagram, group3, {0, 0, 0}, 1, {8, 8}, {}, stroke2);
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "affine frame!!ttt");
    }
  }
  if (contains(tag, {"trianglea"})) {
    auto vertices = dconstants::triangle_positions;
    auto point    = vec3f{0.0, -0.2, 0.0};
    {
      auto group1 = add_group(diagram, {-1.5, 0.0, 0.0});
      add_trianglev(diagram, group1, vertices,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_{2}$!!r",
              "$\\mathbf{p}_{3}$!!r"});
      add_point(diagram, group1, point, "$\\mathbf{p}(w_1,w_2,w_3)$!!b");
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "triangle point!!tt");
    }
    {
      auto group2 = add_group(diagram, {1.5, 0.0, 0.0});
      add_trianglev(diagram, group2, vertices,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_{2}$!!r",
              "$\\mathbf{p}_{3}$!!r"},
          transparent);
      add_point(diagram, group2, point, "$\\mathbf{p}$!!b");
      add_label(
          diagram, group2, {0.0, 1.0, 0.0}, "barycentric coordinates!!tt");
      add_trianglev(diagram, group2, {vertices[0], vertices[1], point},
          {"$w_3$!!e"}, fill2);
      add_trianglev(diagram, group2, {vertices[1], vertices[2], point},
          {"$w_1$!!e"}, fill3);
      add_trianglev(diagram, group2, {vertices[2], vertices[0], point},
          {"$w_2$!!e"}, fill4);
    }
  }
  if (contains(tag, {"line"})) {
    auto vertices = vector<vec3f>{{-0.5, -1.0, 0.0}, {0.5, 1.0, 0.0}};
    auto point    = vec3f{0.0, 0.0, 0.0};
    {
      auto group1 = add_group(diagram, {-1.3, 0.0, 0.0});
      add_line(diagram, group1, vertices[0], vertices[1],
          {"$\\mathbf{p}_1$!!r", "$\\mathbf{p}_2$!!r"});
      add_point(diagram, group1, {0.0, 0.0, 0.0}, "$\\mathbf{p}(u)$!!r");
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "line point!!tt");
    }
    {
      auto group2 = add_group(diagram, {1.4, 0.0, 0.0});
      add_arrow(diagram, group2, vertices[0], vertices[1],
          {"$\\mathbf{p}_1$!!r", "$\\mathbf{e}_{2,1}$!!r"}, black);
      add_point(diagram, group2, point, "$\\mathbf{p}(u)$!!r");
      add_label(
          diagram, group2, {0.0, 1.0, 0.0}, "barycentric coordinates!!tt");
      auto offset = vec3f{-0.2, 0.0, 0.0};
      add_lines(diagram, group2,
          {vertices[0] + offset, vertices[1] + offset, vertices[0],
              vertices[0] + offset, vertices[1] + offset, vertices[1]},
          {}, stroke3);
      add_label(diagram, group2, (vertices[0] + vertices[1]) / 2, "$u$!!l");
    }
  }
  return diagram;
}

diagram_data environment_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto image = load_image(
      "tests/_diagramdata/images/hdrihaven/aerodynamicsworkshop_hdrihaven.hdr");

  if (contains(tag, {"map"})) {
    {
      auto group1 = add_group(diagram, {-1.6, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0},
          "spherical HDR image in lat-long format!!tht");
      add_image(diagram, group1, {0, 0, 0}, 1, image);
      add_arrows(diagram, group1,
          {{-2.1, 1.1, 0.0}, {+2.1, 1.1, 0.0}, {-2.1, 1.1, 0.0},
              {-2.1, -1.1, 0.0}},
          {"", "$s = \\phi / (2\\pi)$!!r", "", "$t = \\theta / \\pi$!!l"});
    }
    {
      auto group2 = add_group(diagram, {2.5, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "environment map!!tht");
      add_uvspheret(
          diagram, group2 * drotation({90.0, 0.0, 215}), {0, 0, 0}, 1, image);
      add_axes(diagram, group2 * drotation({290.0, 0.0, 215}), {0, 0, 0}, 1,
          {"", "$\\mathbf{z}$!!b", "$\\mathbf{x}$!!b", "$\\mathbf{y}$!!l"});
    }
  }
  return diagram;
}

diagram_data shape_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"elements"})) {
    {
      auto group1 = add_group(diagram, {-2.5, 0.0, 0.0});
      add_point(diagram, group1, {0, 0, 0}, "$\\mathbf{p}_1$!!r");
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "points!!tt");
      add_point(diagram, group1, {0.0, 0.0, 0.0}, "", tfill1, extra_thick);
    }
    {
      auto group2 = add_group(diagram, {-0.385, 0.0, 0.0});
      add_line(diagram, group2, {-0.5, -1.0, 0.0}, {+0.5, +1.0, 0.0},
          {"$\\mathbf{p}_1$!!r", "$\\mathbf{p}_2$!!r"});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "lines!!ttt");
      add_line(diagram, group2, {-0.5, -1.0, 0.0}, {+0.5, +1.0, 0.0}, {},
          tfill1, extra_thick * 2);
    }
    {
      auto group3 = add_group(diagram, {2.5, 0.0, 0.0});
      add_triangle(diagram, group3, {0, 0, 0}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "triangles!!tt");
    }
  }
  if (contains(tag, {"elements1"})) {
    auto tetra  = shape_data{.triangles = dconstants::tetra_triangles,
         .positions                     = dconstants::tetra_positions};
    auto labels = vector<string>{"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!b",
        "$\\mathbf{p}_3$!!thl", "$\\mathbf{p}_4$!!r"};
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({20.0, 240.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "points!!tt");
      add_points(diagram, group1, tetra.positions, labels);
      add_points(diagram, group1, tetra.positions, {}, tfill1, extra_thick);
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_shape(diagram, group2 * drotation({20.0, 240.0, 0.0}), tetra, labels,
          black);
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "lines!!ttt");
      add_shape(diagram, group2 * drotation({20.0, 240.0, 0.0}), tetra, {},
          transparent, tfill1, extra_thick * 2);
    }
    {
      auto group3 = add_group(diagram, {3.0, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "triangles!!tt");
      add_shape(diagram, group3 * drotation({20.0, 240.0, 0.0}), tetra, labels,
          tfill1);
    }
  }
  if (contains(tag, {"indexed1"})) {
    {
      auto group1 = add_group(diagram, {-2.0, 0.0, 0.0});
      add_triangles(diagram, group1, {{0, 1, 2}, {1, 2, 3}, {1, 3, 4}},
          {{-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {0.0, 1.0, 0.0},
              {2.0, 1.0, 0.0}, {3.0, -1.0, 0.0}},
          {"$\\mathbf{p}{0}$!!l", "$\\mathbf{p}{1}$!!r", "$\\mathbf{p}{2}$!!l",
              "$\\mathbf{p}{3}$!!r", "$\\mathbf{p}{4}$!!r"});
    }
    {
      auto group2 = add_group(diagram, {2.5, 0.0, 0.0});
      add_labels(diagram, group2,
          {{0.0, 0.4, 0.0}, {0.0, 0.0, 0.0}, {0.0, -0.4, 0.0}, {0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0}},
          {"$\\text{triangle}{0} = (0,1,2)$", "$\\text{triangle}{1} = (1,3,2)$",
              "$\\text{triangle}{2} = (1,4,3)$"});
    }
  }
  if (contains(tag, {"indexed2"})) {
    {
      auto shape  = shape_data{.lines = {{0, 1}, {1, 2}, {2, 3}},
           .positions = {{-1.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, -1.0, 0.0},
               {2.0, 1.0, 0.0}}};
      auto group1 = add_group(diagram, {-1.5, 0.0, 0.0});
      add_lines(diagram, group1, shape,
          {"$\\mathbf{p}{0}$!!l", "$\\mathbf{p}{1}$!!l", "$\\mathbf{p}{2}$!!r",
              "$\\mathbf{p}{3}$!!r"});
      add_lines(diagram, group1, shape, {}, tfill1, extra_thick * 2);
    }
    {
      auto group2 = add_group(diagram, {2.5, 0.0, 0.0});
      add_labels(diagram, group2,
          {{0.0, 0.4, 0.0}, {0.0, 0.0, 0.0}, {0.0, -0.4, 0.0}, {0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0}},
          {"$\\text{line}{0} = (0,1)$", "$\\text{line}{1} = (1,2)$",
              "$\\text{line}{2} = (2,3)$"});
    }
  }
  if (contains(tag, {"indexed3"})) {
    auto tetra  = shape_data{.triangles = dconstants::tetra_triangles,
         .positions                     = dconstants::tetra_positions};
    auto labels = vector<string>{"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!b",
        "$\\mathbf{p}_3$!!thl", "$\\mathbf{p}_4$!!r"};
    {
      auto group1 = add_group(diagram, {-2.0, 0.0, 0.0});
      add_shape(diagram, group1 * drotation({20.0, 240.0, 0.0}), tetra, {},
          transparent);
      add_label(diagram, group1, {0.385, 1.0, 0.0}, "lines!!ttt");
      add_shape(diagram, group1 * drotation({20.0, 240.0, 0.0}), tetra, {},
          transparent, tfill1, extra_thick * 2);
    }
    {
      auto group2 = add_group(diagram, {0.75, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "triangles!!tt");
      add_shape(diagram, group2 * drotation({20.0, 240.0, 0.0}), tetra, labels,
          tfill1);
    }
    {
      auto group3 = add_group(diagram, {-3.5, 0.0, 0.0});
      add_labels(diagram, group3,
          {{0.0, 0.8, 0.0}, {0.0, 0.4, 0.0}, {0.0, 0.0, 0.0}, {0.0, -0.4, 0.0}},
          {"$\\text{line}{0}=(0,1)$", "$\\text{line}{1}=(1,2)$",
              "$\\text{line}{2}=(2,0)$", "$\\cdots$"});
    }
    {
      auto group4 = add_group(diagram, {+3.5, 0.0, 0.0});
      add_labels(diagram, group4,
          {{0.0, 0.8, 0.0}, {0.0, 0.4, 0.0}, {0.0, 0.0, 0.0}, {0.0, -0.4, 0.0}},
          {"$\\text{triangle}{0}=(0,1,2)$", "$\\text{triangle}{1}=(0,1,3)$",
              "$\\text{triangle}{2}=(2,0,3)$",
              "$\\text{triangle}{3}=(1,2,3)$"});
    }
  }
  return diagram;
}

diagram_data shapeapprox_diagrams(const string& tag) {
  auto diagram = make_diagram();

  // Surface approx shapes
  enum struct approx_type { surface, vertex, face };

  auto make_points = [](int steps) {
    auto positions = vector<vec3f>{};
    for (auto idx : range(steps)) {
      auto theta = 2 * pif * idx / (float)steps;
      positions.push_back({yocto::cos(theta), yocto::sin(theta), 0});
    }
    return positions;
  };
  auto dshapeapprox_circle = [&](int steps) {
    auto shape      = shape_data{};
    shape.positions = make_points(steps);
    shape.positions.push_back({0, 0, 0});
    for (auto idx : range(steps))
      shape.triangles.push_back({idx, (idx + 1) % steps, steps});
    for (auto idx : range(steps))
      shape.lines.push_back({idx, (idx + 1) % steps});
    return shape;
  };
  auto dshapeapprox_points = [&](int steps) {
    auto shape      = shape_data{};
    shape.positions = make_points(steps);
    for (auto idx : range(steps)) shape.points.push_back(idx);
    return shape;
  };
  auto dshapeapprox_normals = [&](int steps) {
    auto shape     = shape_data{};
    auto positions = make_points(steps);
    for (auto idx : range(steps)) {
      shape.positions.push_back(positions[idx]);
      shape.positions.push_back(positions[idx] * 1.4f);
    }
    for (auto idx : range(steps))
      shape.lines.push_back({idx * 2 + 0, idx * 2 + 1});
    return shape;
  };
  auto dshapeapprox_fnormals = [&](int steps, int isteps) {
    auto shape     = shape_data{};
    auto positions = make_points(steps);
    for (auto idx : range(steps)) {
      auto p1t = positions[idx], p2t = positions[(idx + 1) % steps];
      auto nt = normalize(p1t - p2t);
      std::swap(nt[0], nt[1]);
      nt[0] = -nt[0];
      if (isteps == 0) {
        shape.positions.push_back((p1t + p2t) / 2);
        shape.positions.push_back((p1t + p2t) / 2 + nt * 0.4f);
      } else {
        for (auto j = 0; j <= isteps; j++) {
          auto p1 = p1t + (p2t - p1t) * (float)j / isteps;
          auto p2 = p1 + nt * 0.4f;
          shape.positions.push_back(p1);
          shape.positions.push_back(p2);
        }
      }
    }
    for (auto idx : range((int)shape.positions.size() / 2))
      shape.lines.push_back({idx * 2 + 0, idx * 2 + 1});
    return shape;
  };
  auto dshapeapprox_vnormals = [&](int steps, int isteps) {
    auto shape     = shape_data{};
    auto positions = make_points(steps);
    for (auto idx : range(steps)) {
      auto p1t = positions[idx], p2t = positions[(idx + 1) % steps];
      auto n1t = normalize(p1t), n2t = normalize(p2t);
      for (auto j = 0; j < isteps; j++) {
        auto p1 = p1t + (p2t - p1t) * (float)j / isteps;
        auto pn = normalize(n1t + (n2t - n1t) * (float)j / isteps);
        auto p2 = p1 + pn * 0.4f;
        shape.positions.push_back(p1);
        shape.positions.push_back(p2);
      }
    }
    for (auto idx : range((int)shape.positions.size() / 2))
      shape.lines.push_back({idx * 2 + 0, idx * 2 + 1});
    return shape;
  };

  auto steps  = 8;
  auto isteps = 4;

  if (contains(tag, {"circle"})) {
    {
      auto group1 = add_group(diagram, {-3.0, 0.0, 0.0});
      add_shape(diagram, group1, dshapeapprox_circle(steps * isteps));
      add_labels(diagram, group1, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
          {"surface!!ht", "normals!!hb"});
      add_arrows(
          diagram, group1, dshapeapprox_normals(steps * isteps).positions);
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_shape(diagram, group2, dshapeapprox_circle(steps));
      add_labels(diagram, group2, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
          {"face!!ht", "normals!!hb"});
      add_arrows(diagram, group2, dshapeapprox_fnormals(steps, 0).positions);
      add_arrows(diagram, group2,
          dshapeapprox_fnormals(steps, isteps).positions, {}, stroke2);
      add_points(
          diagram, group2, dshapeapprox_points(steps).positions, {}, stroke3);
    }
    {
      auto group3 = add_group(diagram, {3.0, 0.0, 0.0});
      add_shape(diagram, group3, dshapeapprox_circle(steps));
      add_labels(diagram, group3, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
          {"vertex!!ht", "normals!!hb"});
      add_arrows(diagram, group3, dshapeapprox_normals(steps).positions);
      add_arrows(diagram, group3,
          dshapeapprox_vnormals(steps, isteps).positions, {}, stroke2);
      add_points(
          diagram, group3, dshapeapprox_points(steps).positions, {}, stroke3);
    }
  }
  return diagram;
}

diagram_data texcoords_diagrams(const string& tag) {
  auto diagram = make_diagram();

  // project a shape onto the uv plane
  auto project_as_texcoords = [](const shape_data& shape_) {
    auto shape = shape_;
    for (auto idx : range(shape.positions.size()))
      shape.positions[idx] = {2 * shape.texcoords[idx].x - 1,
          2 * (1 - shape.texcoords[idx].y) - 1, 0};
    return shape;
  };

  if (contains(tag, {"triangle"})) {
    auto shape   = shape_data{.triangles = {{0, 1, 2}},
          .positions = {{-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {0.0, 1.0, 0.0}},
          .texcoords = {{0.35, 0.4}, {0.65, 0.4}, {0.5, 0.1}}};
    auto uvramp  = load_image("tests/_diagramdata/images/yocto/uvramp.png");
    auto texture = load_image("tests/_diagramdata/images/yocto/uvgrid.png");

    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({15.0, 320.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "texture coordinates!!tt");
      add_shape(diagram, group1, shape,
          {"$(s_1,t_1)$!!b", "$(s_2,t_2)$!!b", "$(s_3,t_3)$!!r"}, uvramp);
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "image texture!!tht");
      add_arrows(diagram, group2,
          {{-1.0, 1.0, 0.0}, {+1.0, +1.0, 0.0}, {-1.0, +1.0, 0.0},
              {-1.0, -1.0, 0.0}},
          {"", "$s$!!b", "", "$t$!!r"});
      add_tquad(diagram, group2, texture);
      add_shape(diagram, group2, project_as_texcoords(shape),
          {"$(s_1,t_1)$!!bl", "$(s_2,t_2)$!!br", "$(s_3,t_3)$!!r"},
          transparent);
    }
    {
      auto group3 = add_group(
          diagram, {3.0, 0.0, 0.0}, drotation({15.0, 320.0, 0.0}));
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "textured triangle!!tht");
      add_shape(diagram, group3, shape,
          {"$(s_1,t_1)$!!b", "$(s_2,t_2)$!!b", "$(s_3,t_3)$!!r"}, texture);
    }
  }
  if (contains(tag, {"cow"})) {
    auto shape = load_shape(
        "tests/_diagramdata/shapes/crane/spot_triangulated.obj");
    auto uvramp  = load_image("tests/_diagramdata/images/yocto/uvramp.png");
    auto texture = load_image(
        "tests/_diagramdata/images/crane/spot_texture.png");

    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({15.0, 130.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "texture coordinates!!tt");
      add_shape(diagram, group1, shape, {}, uvramp);
    }
    {
      auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "image texture!!tht");
      add_arrows(diagram, group2,
          {{-1.0, 1.0, 0.0}, {+1.0, +1.0, 0.0}, {-1.0, +1.0, 0.0},
              {-1.0, -1.0, 0.0}},
          {"", "$s$!!b", "", "$t$!!r"});
      add_tquad(diagram, group2, texture);
      add_shape(diagram, group2, project_as_texcoords(shape), {},
          {0.4, 1.0, 1.0, 0.5});
    }
    {
      auto group3 = add_group(
          diagram, {3.0, 0.0, 0.0}, drotation({15.0, 130.0, 0.0}));
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "textured triangle!!tht");
      add_shape(diagram, group3, shape, {}, texture);
    }
  }
  return diagram;
}

// TODO: remove specials
// Antialiasing shapes
static unordered_map<string, shape_data> builtin_dict() {
  enum struct samples_type { center, supersample, random, stratified };

  auto samples = [](vec2i steps, vec2f size, samples_type type, int ns = 2) {
    auto shape  = shape_data{};
    auto rng    = make_rng(971989);
    auto size3  = vec3f{size.x, size.y, 0};
    auto rand1f = [](auto& rng) { return yocto::rand1f(rng) * 0.8f + 0.1f; };
    for (auto j : range(steps.y)) {
      for (auto i : range(steps.x)) {
        for (auto sj : range(ns)) {
          for (auto si : range(ns)) {
            auto uv = vec2f{0, 0};
            switch (type) {
              case samples_type::center: {
                uv = {(i + 0.5f) / steps.x, (j + 0.5f) / steps.y};
              } break;
              case samples_type::random: {
                uv = {(i + rand1f(rng)) / steps.x, (j + rand1f(rng)) / steps.y};
              } break;
              case samples_type::supersample: {
                uv = {(i + (si + 0.5f) / ns) / steps.x,
                    (j + (sj + 0.5f) / ns) / steps.y};
              } break;
              case samples_type::stratified: {
                uv = {(i + (si + rand1f(rng)) / ns) / steps.x,
                    (j + (sj + rand1f(rng)) / ns) / steps.y};
              } break;
            }
            shape.positions.push_back(
                {size.x * (2 * uv.x - 1), size.y * (2 * uv.y - 1), 0});
            shape.points.push_back((int)shape.positions.size() - 1);
          }
        }
      }
    }
    return shape;
  };
  auto grid = [](vec2i steps, vec2f size) {
    auto shape = make_rect(steps, size);
    return shape_data{
        .lines     = get_edges(shape.quads),
        .positions = shape.positions,
    };
  };
  auto cgrid = [](vec2i steps, vec2f size, const vector<vec3f>& samples,
                   const auto& inside) {
    auto coverage = [](const vector<vec3f>& samples, int samplespp, vec2f size,
                        const auto& inside) {
      auto colors = vector<vec4f>(samples.size() / samplespp, {0, 0, 0, 0});
      for (auto idx : range((int)samples.size())) {
        colors[idx / samplespp] += (inside(samples[idx], size)
                                           ? vec4f{0.4, 1, 1, 1}
                                           : vec4f{1, 1, 1, 1}) /
                                   samplespp;
      }
      return colors;
    };

    auto shape = make_rect(steps, size);
    return shape_data{
        .lines     = get_edges(shape.quads),
        .quads     = shape.quads,
        .positions = shape.positions,
        // .fcolors   = coverage(
        //     samples, (int)samples.size() / (steps.x * steps.y), size,
        //     inside),
    };
  };
  auto poly = [](int steps, const auto& func, vec2f size) {
    auto shape = shape_data{};
    for (auto idx = 0; idx <= steps; idx++) {
      auto u = (2 * idx / (float)steps - 1) * size.x;
      auto v = func(u, size);
      shape.positions.push_back({u, v, 0});
    }
    shape.positions.push_back({-size.x, -size.y, 0});
    for (auto idx : range((int)shape.positions.size() - 1)) {
      shape.triangles.push_back(
          {idx, idx + 1, (int)shape.positions.size() - 1});
    }
    return shape;
  };

  auto func = [](float u, vec2f size) -> float {
    return size.y - size.y * 0.5f * (1 + u / size.x) * (1 + u / size.x);
  };
  auto inside = [](vec3f position, vec2f size) -> bool {
    return position.y <= size.y - size.y * 0.5f * (1 + position.x / size.x) *
                                      (1 + position.x / size.x);
  };

  auto size       = vec2f{1.5, 1};
  auto steps      = vec2i{6, 4};
  auto subsamples = 2;
  auto circle_c   = vec3f{0.15f, 0, 0};
  auto circle_r   = 0.8f;

  auto shapes                           = unordered_map<string, shape_data>{};
  shapes["antialiasing_samples_center"] = samples(
      steps, size, samples_type::center);
  shapes["antialiasing_cgrid_center"]        = cgrid(steps, size,
             samples(steps, size, samples_type::center).positions, inside);
  shapes["antialiasing_samples_supersample"] = samples(
      steps, size, samples_type::supersample);
  shapes["antialiasing_cgrid_supersample"] = cgrid(steps, size,
      samples(steps, size, samples_type::supersample).positions, inside);
  shapes["antialiasing_samples_random"]    = samples(
      steps, size, samples_type::random);
  shapes["antialiasing_cgrid_random"]       = cgrid(steps, size,
            samples(steps, size, samples_type::random).positions, inside);
  shapes["antialiasing_samples_stratified"] = samples(
      steps, size, samples_type::stratified);
  shapes["antialiasing_cgrid_stratified"] = cgrid(steps, size,
      samples(steps, size, samples_type::stratified).positions, inside);
  shapes["antialiasing_grid_pixels"]      = grid(steps, size);
  shapes["antialiasing_grid_supersample"] = grid(steps * 2, size);
  shapes["antialiasing_polygon"]          = poly(32, func, size);
  return shapes;
}

diagram_data antialiasing_diagram(const vector<vec2f>& samples,
    const array2d<vec4f>& image, const array2d<vec4f>& highres) {
  auto ratio = image.extent(0) / (float)image.extent(1);

  auto points = vector<vec3f>{};
  for (auto sample : samples)
    points.push_back({ratio * (sample.x * 2 - 1), sample.y * 2 - 1, 0});

  auto diagram = make_diagram();

  {
    auto group1 = add_group(diagram, {-1.75, 0.0, 0.0});
    add_points(diagram, group1, points);
    add_image(diagram, group1, {0, 0, 0}, 1, highres);
    add_grid(diagram, group1, {0, 0, 0}, 1, (vec2i)image.extents());
    add_label(diagram, group1, {0.0, 1.0, 0.0}, "samples!!t");
  }
  {
    auto group2 = add_group(diagram, {1.75, 0.0, 0.0});
    add_image(diagram, group2, {0, 0, 0}, 1, image);
    add_grid(diagram, group2, {0, 0, 0}, 1, (vec2i)image.extents());
    add_label(diagram, group2, {0.0, 1.0, 0.0}, "reconstruction!!t");
  }

  return diagram;
}

diagram_data antialiasing_diagrams(const string& tag) {
  auto shapes = builtin_dict();

  auto nsamples = 2;
  auto size     = vec2f{1.5, 1};
  auto steps    = vec2i{6, 4};

  auto coverage = [](const vector<vec3f>& samples, int samplespp, vec2f size,
                      const auto& inside) {
    auto colors = vector<vec4f>(samples.size() / samplespp, {0, 0, 0, 0});
    for (auto idx : range((int)samples.size())) {
      colors[idx / samplespp] += (inside(samples[idx], size)
                                         ? vec4f{0.4, 1, 1, 1}
                                         : vec4f{1, 1, 1, 1}) /
                                 samplespp;
    }
    return colors;
  };

  auto shade = [](vec2f sample) {
    auto circle_c = vec3f{0.0f, 1.75f, 0};
    auto circle_r = 1.75f;
    auto point    = vec3f{1.5f * sample.x, sample.y, 0};
    auto inside   = length(point - circle_c) < circle_r;
    return inside ? dcolors::fill1 : vec4f{1.0, 1.0, 1.0, 1.0};
  };

  auto rand1f = [](auto& rng) { return yocto::rand1f(rng) * 0.8f + 0.1f; };

  auto highres = array2d<vec4f>({640, 480});
  for (auto ij : range(highres.extents())) {
    auto uv     = (ij + 0.5f) / highres.extents();
    highres[ij] = shade(uv);
  }

  auto rng     = make_rng(971989);
  auto samples = vector<vec2f>{};
  auto image   = array2d<vec4f>({6, 4});

  if (contains(tag, {"center"})) {
    for (auto ij : range(steps)) {
      auto uv = (ij + 0.5f) / highres.extents();
      samples.push_back(uv);
      image[ij] = shade(uv);
    }
  }
  if (contains(tag, {"random"})) {
    for (auto ij : range(steps)) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto uv = (ij + rand2f(rng)) / steps;
        samples.push_back(uv);
        image[ij] += shade(uv) / 4;
      }
    }
  }
  if (contains(tag, {"stratified"})) {
    for (auto ij : range(steps)) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto uv = (ij + (sij + rand2f(rng)) / nsamples) / steps;
        samples.push_back(uv);
        image[ij] += shade(uv) / 4;
      }
    }
  }
  if (contains(tag, {"supersample"})) {
    for (auto ij : range(steps)) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto uv = (ij + (sij + 0.5f) / nsamples) / steps;
        samples.push_back(uv);
        image[ij] += shade(uv) / 4;
      }
    }
  }

  return antialiasing_diagram(samples, image, highres);
}

diagram_data brdfframe_diagrams(const string& tag) {
  auto diagram = make_diagram();

  {
    auto group = add_group(
        diagram, {0, 0, 0}, drotation({15.0, 40.0, 0.0}) * yup3x4f);
    {
      auto surface = group;
      if (!contains(tag, {"emission"})) {
        add_quad(diagram, surface, {0, 0, 0}, 0.75, {}, tfill1);
      }
      if (contains(tag, {"emission"})) {
        add_quad(diagram, surface, {0, 0, 0}, 0.75, {}, tfill5);
      }
      add_vector(diagram, surface, {0.0, 0.0, 1.0}, {"", "$\\mathbf{n}$!!bhl"});
      add_point(diagram, surface, {0.0, 0.0, 0.0}, "$\\mathbf{p}$!!lln");
    }
    {
      auto outgoing = group;
      add_vector(diagram, outgoing, {-0.7, 0.0, 0.7}, {"", "$\\mathbf{o}$!!b"});
      add_line(
          diagram, outgoing, {-1.25, 0.0, 1.25}, {-0.8, 0.0, 0.8}, {}, gray);
      add_disk(diagram, outgoing * dlookat({-1.25, 0.0, 1.25}), {0, 0, 0}, 0.2,
          {}, tfill3);
      add_point(diagram, outgoing, {-1.25, 0.0, 1.25}, "$\\mathbf{e}$!!ll");
    }
    {
      auto incoming = group;
      if (contains(tag, {"frame"})) {
        auto frame = incoming;
        add_vector(diagram, frame, {0.7, 0.0, 0.7}, {"", "$\\mathbf{i}$!!b"});
        add_line(diagram, frame, {1.25, 0.0, 1.25}, {0.8, 0.0, 0.8}, {}, gray);
      }
      if (contains(tag, {"arealight", "areashadow", "diffuse", "diffuse2",
                            "envlight", "envshadow"})) {
        auto random = incoming;
        add_vector(diagram, random, {0.7, 0.0, 0.7},
            {"", "$\\mathbf{i} = \\text{random}(\\mathbf{n})$!!bar"});
        if (!contains(tag, {"areashadow", "envlight", "envshadow"})) {
          add_line(
              diagram, random, {1.25, 0.0, 1.25}, {0.8, 0.0, 0.8}, {}, gray);
        }
        if (contains(tag, {"envlight"})) {
          add_line(diagram, random, {1.5, 0.0, 1.5}, {0.8, 0.0, 0.8}, {}, gray);
        }
        if (contains(tag, {"areashadow", "envshadow"})) {
          add_line(diagram, random, {0.9, 0.0, 0.9}, {0.8, 0.0, 0.8}, {}, gray);
        }
        add_rlines(diagram, random, {0, 0, 0}, 1, rlines_type::hemicos, 8,
            false, {}, gray);
      }
      if (contains(tag, {"areashadow2", "pointshadow"})) {
        auto occluded = incoming;
        add_vector(diagram, occluded, {0.7, 0.0, 0.7},
            {"",
                "$\\mathbf{i} = \\text{normalize}(\\mathbf{s} - \\mathbf{p})$!!bar"});
        add_line(
            diagram, occluded, {0.95, 0.0, 0.95}, {0.8, 0.0, 0.8}, {}, gray);
      }
      if (contains(tag, {"arealight2", "pointlight"})) {
        auto light = incoming;
        add_vector(diagram, light, {0.7, 0.0, 0.7},
            {"",
                "$\\mathbf{i} = \\text{normalize}(\\mathbf{s} - \\mathbf{p})$!!bar"});
        add_line(diagram, light, {1.25, 0.0, 1.25}, {0.8, 0.0, 0.8}, {}, gray);
      }
      if (contains(tag, {"areashadow2", "pointshadow"})) {
        auto light_occluded = incoming;
        add_vector(diagram, light_occluded, {0.7, 0.0, 0.7},
            {"",
                "$\\mathbf{i} = \\text{normalize}(\\mathbf{s} - \\mathbf{p})$!!bar"});
      }
      if (contains(tag, {"reflection", "transmission"})) {
        auto reflection = incoming;
        add_vector(diagram, reflection, {0.7, 0.0, 0.7},
            {"",
                "$\\mathbf{i} = \\text{reflect}(\\mathbf{o},\\mathbf{n}) = - \\mathbf{o} + 2(\\mathbf{o} \\cdot \\mathbf{n}) \\mathbf{o}$!!bar"});
        add_line(
            diagram, reflection, {1.25, 0.0, 1.25}, {0.8, 0.0, 0.8}, {}, gray);
      }
      if (contains(tag, {"rreflection"})) {
        auto rreflection = incoming;
        add_vector(diagram, rreflection, {0.7, 0.0, 0.7},
            {"",
                "$\\mathbf{h} = \\text{random}(\\mathbf{n},\\alpha) \\ with \\  \\mathbf{i} = \\text{reflect}(\\mathbf{o}, \\mathbf{h})$!!bar"});
        add_line(
            diagram, rreflection, {1.25, 0.0, 1.25}, {0.8, 0.0, 0.8}, {}, gray);
        add_rlines(diagram,
            rreflection * drotation({45.0, 69.23100280761719, 0.0}), {0, 0, 0},
            1, rlines_type::hemicospower, 4, false, {}, gray);
      }
      if (contains(tag, {"transmission"})) {
        auto transmission = incoming;
        add_vector(diagram, transmission, {0.7, 0.0, -0.7},
            {"",
                "$\\mathbf{i} = \\text{transmit}(\\mathbf{o},\\mathbf{n}) = - \\mathbf{o}$!!tar"});
        add_line(diagram, transmission, {1.25, 0.0, -1.25}, {0.8, 0.0, -0.8},
            {}, gray);
      }
      if (contains(tag, {"diffuse2"})) {
        auto random2 = incoming;
        add_vector(diagram, random2, {0.7, 0.7, 0.2},
            {"", "$\\mathbf{i}=\\cdots$!!hbr"});
        add_line(diagram, random2, {0.8, 0.8, 0.2}, {1.5, 1.5, 0.5}, {}, gray);
      }
    }
    if (contains(tag, {"diffuse", "diffuse2", "reflection", "rreflection",
                          "transmission", "arealight", "areashadow"})) {
      auto next = group * dlookat({1.25, 0.0, 1.25});
      if (!contains(tag, {"arealight", "areashadow"})) {
        add_quad(diagram, next, {0, 0, 0}, 0.4, {}, tfill2);
      }
      if (!contains(tag, {"areashadow"})) {
        add_point(diagram, next, {0.0, 0.0, 0.0},
            "$\\mathbf{q} = \\text{raytrace}(\\mathbf{p}, \\mathbf{i})$!!rrr");
      }
      if (contains(tag, {"arealight", "areashadow"})) {
        add_quad(diagram, next, {0, 0, 0}, 0.4, {}, tfill5);
      }
    }
    if (contains(tag, {"arealight2", "areashadow2"})) {
      auto arealight = group * dlookat({1.25, 0.0, 1.25});
      add_quad(diagram, arealight, {0, 0, 0}, 0.4, {}, tfill5);
      add_point(diagram, arealight, {0.0, 0.0, 0.0},
          "$\\mathbf{s} = \\text{random}(surface)$!!rrr");
      if (contains(tag, {"arealight2", "areashadow2"})) {
        add_rpoints(diagram, arealight, {0, 0, 0}, 0.4, 16, true, {}, gray);
      }
    }
    if (contains(tag, {"frame", "pointlight", "pointshadow"})) {
      auto pointlight = group * dlookat({1.25, 0.0, 1.25});
      add_disk(diagram, pointlight, {0, 0, 0}, 0.2, {}, tfill5);
      add_point(diagram, pointlight, {0.0, 0.0, 0.0}, "$\\mathbf{s}$!!rr");
    }
    if (contains(tag, {"envlight", "envshadow"})) {
      auto envlight = group;
      add_sphere(diagram, envlight, {0, 0, 0}, 2, {}, {1.0, 1.0, 0.0, 0.2});
      add_point(diagram, envlight, {1.5, 0.0, 1.5},
          "$\\text{environment}(\\mathbf{i})$!!rr");
    }
    if (contains(
            tag, {"areashadow", "areashadow2", "pointshadow", "envshadow"})) {
      auto occluder = group * dlookat({0.95, 0.0, 0.95});
      add_quad(diagram, occluder, {0, 0, 0}, 0.4, {}, tfill2);
      add_point(diagram, occluder, {0.0, 0.0, 0.0},
          "$\\mathbf{q} = \\text{raytrace}(\\mathbf{p},\\mathbf{i})$!!rrr");
    }
    if (contains(tag, {"diffuse2"})) {
      auto next2 = group * dlookat({1.5, 1.5, 0.5});
      add_quad(diagram, next2, {0.0, 0.0, 0.0}, 0.4, {}, tfill2);
      add_point(diagram, next2, {0.0, 0.0, 0.0},
          "$\\mathbf{q} = \\text{raytrace}(\\mathbf{p}, \\mathbf{i})$!!rr");
    }
    if (contains(tag, {"transmission"})) {
      auto nextt = group * dlookat({1.25, 0.0, -1.25});
      add_quad(diagram, nextt, {0, 0, 0}, 0.4, {}, tfill2);
      add_point(diagram, nextt, {0.0, 0.0, 0.0},
          "$\\mathbf{q} = \\text{raytrace}(\\mathbf{p}, \\mathbf{i})$!!rrr");
    }
  }
  return diagram;
}

auto plot_ggx(float angle, float roughness) -> float {
  auto cosine     = cos(angle);
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  return roughness2 / (pif * (cosine2 * roughness2 + 1 - cosine2) *
                          (cosine2 * roughness2 + 1 - cosine2));
}
auto plot_phong_d(float angle, float roughness) -> float {
  auto cosine   = clamp(cos(angle), 0.0001f, 1.0f);
  auto exponent = 2 / (roughness * roughness) - 2;
  return (exponent + 2) / (2 * pif) * pow(cosine, exponent);
}
auto plot_g1(float angle, float roughness) -> float {
  auto cosine     = cos(angle);
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  return 2 * abs(cosine) /
         (abs(cosine) + sqrt(cosine2 - roughness2 * cosine2 + roughness2));
};
auto plot_fresnels(float angle, float specular) -> float {
  auto cosine = cos(angle);
  return specular +
         (1 - specular) * pow(clamp(1 - abs(cosine), 0.0f, 1.0f), 5.0f);
}
auto plot_fresneld(float angle, float eta) -> float {
  auto cosw = cos(angle);

  auto sin2 = 1 - cosw * cosw;
  auto eta2 = eta * eta;

  auto cos2t = 1 - sin2 / eta2;
  if (cos2t < 0) return 1.0f;  // tir

  auto t0 = sqrt(cos2t);
  auto t1 = eta * t0;
  auto t2 = eta * cosw;

  auto rs = (cosw - t1) / (cosw + t1);
  auto rp = (t0 - t2) / (t0 + t2);

  return (rs * rs + rp * rp) / 2;
}
auto plot_ggx_polar(float wo, float wi, float r, bool all = true) -> float {
  auto wh     = (wo + wi) / 2;
  auto ndo    = cos(wo);
  auto ndi    = cos(wi);
  auto ndh    = cos(wh);
  auto cos2   = ndh * ndh;
  auto tan2   = (1 - cos2) / cos2;
  auto alpha2 = r * r;
  auto d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
  auto lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
  auto lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
  auto g        = 1 / (1 + lambda_o + lambda_i);
  if (all)
    return max(0.0f, (d * g / (4 * ndi * ndo)));
  else
    return max(0.0f, (d * g / (4 * ndo)));
}

diagram_data brdfplot_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"fdielectric"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"diectric Fresnel"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_fresneld(x, 1.5); },
        {0.0, pif / 2}, {}, stroke1);
  }
  if (contains(tag, {"ggx"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"ggx microfacet distribution"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ybounds = {0.0, 8.0},
             .ylabels = {"0", "8"},
             .yticks  = {0.0, 8.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_ggx(x, 0.2); }, {0.0, pif / 2},
        {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return plot_ggx(x, 0.3); }, {0.0, pif / 2},
        {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return plot_ggx(x, 0.4); }, {0.0, pif / 2},
        {}, stroke3);
  }
  if (contains(tag, {"phong"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"Phong microfacet distribution"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ybounds = {0.0, 8.0},
             .ylabels = {"0", "8"},
             .yticks  = {0.0, 8.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_phong_d(x, 0.2); },
        {0.0, pif / 2}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return plot_phong_d(x, 0.3); },
        {0.0, pif / 2}, {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return plot_phong_d(x, 0.4); },
        {0.0, pif / 2}, {}, stroke3);
  }
  if (contains(tag, {"schlickd"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"Schlick Fresnel for dielectics"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_fresnels(x, 0.04); },
        {0.0, pif / 2}, {}, stroke1);
  }
  if (contains(tag, {"schlickm"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"Schlick Fresnel for dielectics"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_fresnels(x, 0.4); },
        {0.0, pif / 2}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return plot_fresnels(x, 0.5); },
        {0.0, pif / 2}, {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return plot_fresnels(x, 0.8); },
        {0.0, pif / 2}, {}, stroke3);
  }
  if (contains(tag, {"shadowing"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.aspect     = {pif / 2, 1.0},
             .titles  = {"microfacet shadowing"},
             .xbounds = {0.0, pif / 2},
             .xlabels = {"0", "$\\pi/2$"},
             .xticks  = {0.0, pif / 2},
             .ybounds = {0.0, 1.1},
             .ylabels = {"0", "1"},
             .yticks  = {0.0, 1.0}});
    add_lineplot(
        diagram, plot, [](float x) { return plot_g1(x, 0.01); }, {0.0, pif / 2},
        {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return plot_g1(x, 0.1); }, {0.0, pif / 2},
        {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return plot_g1(x, 0.25); }, {0.0, pif / 2},
        {}, stroke3);
  }
  return diagram;
}

diagram_data cameraray_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"pinhole"})) {
    auto group = add_group(diagram, {0, 0, 0}, drotation({20.0, 300.0, 0.0}));
    add_axes(diagram, group, {0, 0, 0}, 1,
        {"$\\mathbf{o}$!!bl", "$\\mathbf{x}$!!l", "$\\mathbf{y}$!!l",
            "$\\mathbf{z}$!!t"});
    add_disk(diagram, group, {0, 0, 0}, 0.2, {}, tfill3);
    add_rect(diagram, group, {3, 2}, {0.0, 0.0, 1.5}, 0.4,
        {"", "", "", "film!!l"}, tfill1);
    add_rect(diagram, group, {3, 2}, {0.0, 0.0, -3.0}, 0.8,
        {"", "focus plane!!hbhl"}, tfill2);
    add_cline(diagram, group, {0.0, 0.0, 1.1}, {0.0, 0.0, 1.5}, {}, gray);
    add_point(diagram, group, {0.3, -0.3, 1.5}, "$\\mathbf{q}(u,v)$!!bl");
    add_cline(diagram, group, {0.3, -0.3, 1.5}, {0.0, 0.0, 0.0}, {}, gray);
    add_arrow(diagram, group, {0, 0, 0}, {-0.16, 0.16, -0.8},
        {"$\\mathbf{e}$!!r", "$\\mathbf{d}$!!b"}, black);
    add_cline(diagram, group, {-0.16, 0.16, -0.8}, {-0.8, 0.8, -4.0}, {}, gray);
  }
  return diagram;
}

diagram_data lighting_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"lambertlaw"})) {
    {
      auto group1 = add_group(
          diagram, {-1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}) * yup3x4f);
      add_disk(diagram, group1 * dlookat({0.0, 0.0, 1.15}), {0, 0, 0}, 0.5, {},
          tfill5);
      add_rlines(diagram,
          group1 * dlookat({0.0, 0.0, 1.15}) * dscaling({1.0, 1.0, 2.5}),
          {0, 0, 0}, 1, rlines_type::beam, 16, true, {}, stroke5);
      add_disk(diagram, group1, {0, 0, 0}, 0.5, {}, transparent);
      add_quad(diagram, group1, {0, 0, 0}, 1);
    }
    {
      auto group2 = add_group(
          diagram, {1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}) * yup3x4f);
      add_disk(diagram, group2 * dlookat({0.0, -1.0, 1.0}), {0, 0, 0}, 0.5, {},
          tfill5);
      add_rlines(diagram,
          group2 * dlookat({0.0, -1.0, 1.0}) * dscaling({1.0, 1.0, 3.33}),
          {0, 0, 0}, 1, rlines_type::beam, 16, true, {}, stroke5);
      add_disk(diagram, group2 * dscaling({1.0, 1.4199999570846558, 1.0}),
          {0, 0, 0}, 0.5, {}, transparent);
      add_quad(diagram, group2, {0, 0, 0}, 1);
    }
  }
  if (contains(tag, {"lambertlaw2"})) {
    {
      auto group1 = add_group(
          diagram, {-1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}) * yup3x4f);
      add_disk(diagram, group1 * dlookat({0.0, 0.0, 1.15}), {0, 0, 0}, 0.5, {},
          tfill5);
      add_rlines(diagram,
          group1 * dlookat({0.0, 0.0, 1.15}) * dscaling({1.0, 1.0, 2.5}),
          {0, 0, 0}, 1, rlines_type::beam, 16, true, {}, stroke5);
      add_disk(diagram, group1, {0, 0, 0}, 0.5, {}, transparent);
      add_quad(diagram, group1, {0, 0, 0}, 1);
    }
    {
      auto group2 = add_group(
          diagram, {1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}) * yup3x4f);
      add_disk(diagram, group2 * dlookat({0.0, -1.0, 1.0}), {0, 0, 0}, 0.5, {},
          tfill5);
      add_rlines(diagram,
          group2 * dlookat({0.0, -1.0, 1.0}) * dscaling({1.0, 1.0, 3.33}),
          {0, 0, 0}, 1, rlines_type::beam, 16, true, {}, stroke5);
      add_disk(diagram, group2 * dscaling({1.0, 1.4199999570846558, 1.0}),
          {0, 0, 0}, 0.5, {}, transparent);
      add_quad(diagram, group2, {0, 0, 0}, 1);
    }
  }
  if (contains(tag, {"lightrsquared"})) {
    {
      auto group1 = add_group(
          diagram, {-1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}));
      add_point(diagram, group1, {0.0, 0.0, 0.0});
      add_line(diagram, group1, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
          {"r=1!!hthrc"}, gray);
      add_sphere(diagram, group1, {0, 0, 0}, 1, {}, {1.0, 1.0, 0.0, 0.3});
      add_point(diagram, group1, {1.0, 0.0, 0.0}, "$I$!!r");
    }
    {
      auto group2 = add_group(
          diagram, {1.7, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}));
      add_point(diagram, group2, {0.0, 0.0, 0.0});
      add_line(diagram, group2, {0.0, 0.0, 0.0}, {1.5, 0.0, 0.0}, {"r!!hthrc"},
          gray);
      add_sphere(diagram, group2, {0, 0, 0}, 1.5, {}, {1.0, 1.0, 0.0, 0.2});
      add_point(diagram, group2, {1.5, 0.0, 0.0}, "$I/r^2$!!r");
    }
  }
  return diagram;
}

diagram_data rendering_diagrams(const string& tag) {
  auto diagram = make_diagram();

  {
    auto group = add_group(diagram, {0, 0, 0}, drotation({20.0, 300.0, 0.0}));
    add_disk(diagram, group, {0, 0, 0}, 0.2, {}, tfill3);
    add_label(diagram, group, {0.0, 0.0, 0.0}, "lens!!tlhl");
    add_rect(diagram, group, {3, 2}, {0.0, 0.0, 1.5}, 0.4,
        {"", "", "", "film!!l"}, tfill1);
    add_point(diagram, group, {0.3, -0.3, 1.5}, "pixel!!bl");
    add_cube(diagram, group, {-0.2, 0.1, -3.0}, 0.5, {}, fill2);
    {
      auto hit1 = group * dtranslation({-0.5, 0.5, -2.5});
      add_point(diagram, hit1, {0.0, 0.0, 0.0}, "visible \\\\ point!!tt");
      add_sphere(diagram, hit1, {0, 0, 0}, 0.1, {}, {1.0, 1.0, 0.0, 0.5});
    }
    add_cube(diagram, group, {-1.5, 1.1, -0.74}, 0.5, {}, fill2);
    {
      auto hit2 = group * dtranslation({-1.0, 0.75, -1.0});
      add_point(diagram, hit2, {0.0, 0.0, 0.0}, "visible \\\\ point!!bb");
      add_sphere(diagram, hit2, {0, 0, 0}, 0.1, {}, {1.0, 1.0, 0.0, 0.5});
    }
    {
      auto light = group * dtransform({-0.5, 1.5, -1.5}, {0.0, 160, 73});
      add_disk(diagram, light * yup3x4f, {0, 0, 0}, 0.1, {}, fill5);
      add_label(diagram, light, {0.0, 0.0, 0.0}, "light!!tr");
    }
    if (contains(tag, {"image"})) {
      add_grid(diagram, group, {0.0, 0.0, 1.49}, 0.4, {6, 4}, {}, gray);
    }
    if (contains(tag, {"image"})) {
      auto forward_rays = group;
      add_carrow(diagram, forward_rays, {0.3, -0.3, 1.5}, {0.0, 0.0, 0.0}, {},
          stroke4);
      add_carrow(diagram, forward_rays, {0.0, 0.0, 0.0}, {-0.5, 0.5, -2.5}, {},
          stroke4);
      add_carrow(diagram, forward_rays, {-0.5, 0.5, -2.5}, {-0.5, 1.5, -1.5},
          {}, stroke4);
      add_carrow(diagram, forward_rays, {-0.5, 0.5, -2.5}, {-1.0, 0.75, -1.0},
          {}, stroke4);
      add_carrow(diagram, forward_rays, {-1.0, 0.75, -1.0}, {-0.5, 1.5, -1.5},
          {}, stroke4);
    }
    if (contains(tag, {"lighting"})) {
      auto backward_rays = group;
      add_carrow(diagram, backward_rays, {0.0, 0.0, 0.0}, {0.3, -0.3, 1.5}, {},
          stroke4);
      add_carrow(diagram, backward_rays, {-0.5, 0.5, -2.5}, {0.0, 0.0, 0.0}, {},
          stroke4);
      add_carrow(diagram, backward_rays, {-0.5, 1.5, -1.5}, {-0.5, 0.5, -2.5},
          {}, stroke4);
      add_carrow(diagram, backward_rays, {-1.0, 0.75, -1.0}, {-0.5, 0.5, -2.5},
          {}, stroke4);
      add_carrow(diagram, backward_rays, {-0.5, 1.5, -1.5}, {-1.0, 0.75, -1.0},
          {}, stroke4);
    }
    if (contains(tag, {"lighting"})) {
      auto sparkles = group;
      add_rlines(diagram,
          sparkles * dtransform({-0.5, 0.5, -2.5}, {0.0, 0.0, 100.0}),
          {0, 0, 0}, 0.25, rlines_type::hemi, 16, true, {}, stroke5);
      add_rlines(diagram,
          sparkles * dtransform({-1.0, 0.75, -1.0}, {0.0, 0.0, 270.0}),
          {0, 0, 0}, 0.25, rlines_type::hemi, 16, true, {}, stroke5);
      add_rlines(diagram,
          sparkles * dtransform({-0.5, 1.5, -1.5}, {0.0, 160, 75}), {0, 0, 0},
          0.25, rlines_type::hemi, 16, true, {}, stroke5);
    }
  }
  return diagram;
}

diagram_data bbox_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"slabs"})) {
    {
      auto group1 = add_group(
          diagram, {-2.75, 0.0, 0.0}, drotation({15.0, 295.0, 0.0}));
      add_cube(diagram, group1, {0, 0, 0}, 0.49);
      add_slabs(diagram, group1, {0, 0, 0}, 0.5, {}, tfill2, dcolors::stroke2);
    }
  }
  return diagram;
}

diagram_data intersect_diagrams(const string& tag) {
  auto diagram = make_diagram();

  {
    auto group = add_group(diagram, {0, 0, 0}, drotation({20.0, 300.0, 0.0}));
    if (!contains(tag, {"instance"})) {
      add_ray(diagram, group, {0.0, 0.0, 2.5}, {0.0, 0.0, -1.0},
          {"$\\mathbf{e}$!!b", "$\\mathbf{d}$!!b"});
    }
    if (!contains(tag, {"instance"})) {
      add_cline(diagram, group, {0.0, 0.0, 3.5}, {0.0, 0.0, 6.0}, {}, gray);
    }
    if (contains(tag, {"ray"})) {
      add_line(diagram, group, {0.0, 0.0, 0.0}, {0.0, 0.0, 2.5}, {"$t$!!tc"},
          stroke3);
    }
    if (contains(tag, {"triangle"})) {
      auto triangle = group;
      add_triangle(diagram, triangle, {0, 0, 0}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"});
      add_point(diagram, triangle, {0.0, 0.0, 0.0}, "$\\mathbf{r}(t)$!!tnhn");
      add_point(diagram, triangle, {0.0, 0.0, 0.0}, "$\\mathbf{p}(u,v)$!!bnhn");
    }
    if (contains(tag, {"line"})) {
      auto line = group;
      add_line(diagram, line, {0.2, -1.0, 0.0}, {0.2, 1.0, 0.0},
          {"$\\mathbf{p}_1$!!r", "$\\mathbf{p}_2$!!r"});
      if (contains(tag, {"line"})) {
        add_line(diagram, line, {0.2, -1.0, 0.0}, {0.2, 1.0, 0.0}, {}, tfill1,
            extra_thick);
      }
      add_points(diagram, line, {{0.0, 0.0, 0.0}, {0.2, 0.0, 0.0}},
          {"$\\mathbf{r}(t)$!!tl", "$\\mathbf{p}(s)$!!r"});
    }
    if (contains(tag, {"point"})) {
      auto point = group;
      add_point(diagram, point, {0.2, 0.0, 0.0}, "$\\mathbf{p}_1$!!r");
      add_point(diagram, point, {0.2, 0.0, 0.0}, "", tfill1, extra_thick);
      add_point(diagram, point, {0.0, 0.0, 0.0}, "$\\mathbf{r}(t)$!!tl");
    }
    if (contains(tag, {"bbox"})) {
      auto box = group * dtranslation({0.0, 0.0, -1.5});
      add_cube(diagram, box, {0, 0, 0}, 1, {}, tfill1);
      add_labels(diagram, box, {{-1.0, -1.0, -1.5}, {+1.0, +1.0, +1.0}},
          {"$(x_m,y_m,z_m)$!!bhl", "$(x_M,y_M,z_M)$!!thr"});
      add_point(diagram, box, {0.0, 0.0, +1.0});
    }
    if (contains(tag, {"slab"})) {
      auto slab = group;
      add_slab(diagram, slab, {0, 0, 0}, 0.5,
          {"", "$z_1$!!r", "", "", "$z_2$!!r"}, tfill1, dcolors::stroke1);
      add_points(diagram, slab, {{0.0, 0.0, 0.5}, {0.0, 0.0, -0.5}},
          {"$\\mathbf{r}(t_1)$!!tnhn", "$\\mathbf{r}(t_2)$!!tnhn"});
    }
    if (contains(tag, {"shape"})) {
      auto shape = group;
      add_triangle(diagram, shape, {0, 0, 0}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"});
      add_point(diagram, shape, {0.0, 0.0, 0.0}, "$\\mathbf{r}(t_1)$!!bnhn");
      add_triangle(diagram, shape, {0.0, 0.0, -1.852}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"},
          fill2);
      add_point(diagram, shape, {0.0, 0.0, -1.852}, "$\\mathbf{r}(t_2)$!!bnhn");
    }
    if (contains(tag, {"scene"})) {
      auto scene_ = group;
      add_triangle(diagram, scene_, {0, 0, 0}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"});
      add_point(diagram, scene_, {0.0, 0.0, 0.0}, "$\\mathbf{r}(t_1)$!!bnhn");
      add_triangle(diagram, scene_, {0.0, 0.0, -1.852}, 1,
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"},
          fill2);
      add_point(
          diagram, scene_, {0.0, 0.0, -1.852}, "$\\mathbf{r}(t_2)$!!bnhn");
    }
  }
  {
    if (contains(tag, {"instance"})) {
      {
        auto group1 = add_group(
            diagram, {-1.5, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}));
        add_ray(diagram, group1, {0.0, 0.0, 2.5}, {0.0, 0.0, -1.0},
            {"$\\mathbf{e}$!!b", "$\\mathbf{d}$!!b"});
        add_triangle(diagram, group1, {0, 0, 0}, 1,
            {"$\\mathbf{F} \\cdot \\mathbf{p}_1$!!l",
                "$\\mathbf{F} \\cdot \\mathbf{p}_2$!!r",
                "$\\mathbf{F} \\cdot \\mathbf{p}_3$!!r"});
        add_cline(
            diagram, group1, {0.0, 0.0, -1.0}, {0.0, 0.0, -4.0}, {}, gray);
        add_point(diagram, group1, {0.0, 0.0, 0.0}, "$\\mathbf{r}(t)$!!bnhn");
      }
      {
        auto group2 = add_group(
            diagram, {2.5, 0.0, 0.0}, drotation({20.0, 300.0, 0.0}));
        add_ray(diagram, group2, {0.0, 0.0, 2.5}, {0.0, 0.0, -1.0},
            {"$\\mathbf{e}$!!b", "$\\mathbf{d}$!!b"});
        add_triangle(diagram, group2, {0, 0, 0}, 1,
            {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!r"});
        add_cline(
            diagram, group2, {0.0, 0.0, -1.0}, {0.0, 0.0, -4.0}, {}, gray);
        add_point(diagram, group2, {0.0, 0.0, 0.0},
            "$\\mathbf{F}^{-1} \\cdot \\mathbf{r}(t)$!!rhrnn");
      }
    }
  }
  return diagram;
}

diagram_data bezier_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"splines", "hulls", "splits", "tangents"})) {
    {
      auto group1 = add_group(diagram, {-1.5, 0.0, 0.0});
      add_qbezier(diagram, group1,
          {{-1.0, -0.75, 0.0}, {0.0, 0.75, 0.0}, {1.0, 0.0, 0.0}});
      add_polyline(diagram, group1,
          {{-1.0, -0.75, 0.0}, {0.0, 0.75, 0.0}, {1.0, 0.0, 0.0}},
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!l", "$\\mathbf{p}_3$!!r"},
          stroke1);
      add_point(diagram, group1, {0.0, 0.1875, 0.0}, "$\\mathbf{p}(u)$!!b");
      if (contains(tag, {"hulls"})) {
        add_trianglev(diagram, group1,
            {{-1.0, -0.75, 0.0}, {0.0, 0.75, 0.0}, {1.0, 0.0, 0.0}}, {}, fill3);
      }
      if (contains(tag, {"tangents"})) {
        add_arrow(diagram, group1, {-1.0, -0.75, 0.0}, {-0.45, 0.08, 0.0},
            {"", "$\\mathbf{p}_2 - \\mathbf{p}_1$!!l"}, stroke2, 0.0155);
        add_arrow(diagram, group1, {1.0, 0.0, 0.0}, {0.2, 0.6, 0.0},
            {"", "$\\mathbf{p}_2 - \\mathbf{p}_3$!!r"}, stroke2, 0.0155);
      }
      if (contains(tag, {"splits"})) {
        add_polyline(diagram, group1, {{-0.5, 0.0, 0.0}, {0.5, 0.375, 0.0}},
            {"$\\mathbf{q}_1(t)$!!l", "$\\mathbf{q}_2(t)$!!r"}, stroke2);
      }
    }
    {
      auto group2 = add_group(diagram, {1.5, 0.0, 0.0});
      add_bezier(diagram, group2,
          {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
              {1.0, 0.0, 0.0}});
      add_point(diagram, group2, {-0.0375, 0.43, 0.0}, "$\\mathbf{p}(u)$!!b");
      add_polyline(diagram, group2,
          {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
              {1.0, 0.0, 0.0}},
          {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!l", "$\\mathbf{p}_3$!!r",
              "$\\mathbf{p}_4$!!r"},
          stroke1);
      if (contains(tag, {"hulls"})) {
        add_quadv(diagram, group2,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}},
            {}, fill3);
      }
      if (contains(tag, {"tangents"})) {
        add_arrow(diagram, group2, {-1.0, -0.75, 0.0}, {-0.58, 0.16, 0.0},
            {"", "$\\mathbf{p}_2 - \\mathbf{p}_1$!!l"}, stroke3);
        add_arrow(diagram, group2, {1.0, 0.0, 0.0}, {0.29, 0.71, 0.0},
            {"", "$\\mathbf{p}_3 - \\mathbf{p}_4$!!r"}, stroke3);
      }
      if (contains(tag, {"splits"})) {
        add_polyline(diagram, group2,
            {{-0.68, -0.05, 0.0}, {-0.05, 0.7, 0.0}, {0.625, 0.375, 0.0}},
            {"$\\mathbf{q}_1(t)$!!l", "$\\mathbf{q}_2(t)$!!t",
                "$\\mathbf{q}_3(t)$!!tr"},
            stroke2);
      }
      if (contains(tag, {"splits"})) {
        add_line(diagram, group2, {-0.36, 0.32, 0.0}, {0.29, 0.54, 0.0},
            {"$\\mathbf{r}_1(t)$!!l", "$\\mathbf{r}_2(t)$!!r"}, stroke3);
      }
    }
  }
  {
    if (contains(tag, {"cusps"})) {
      {
        auto group1 = add_group(diagram, {-2.0, 0.0, 0.0});
        add_bezier(diagram, group1,
            {{-1.0, 0.0, 0.0}, {-0.25, -0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}});
        add_polyline(diagram, group1,
            {{-1.0, 0.0, 0.0}, {-0.25, -0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}},
            {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!l", "$\\mathbf{p}_3$!!r",
                "$\\mathbf{p}_4$!!r"},
            stroke1);
      }
      {
        auto group2 = add_group(diagram, {2.0, 0.0, 0.0});
        add_bezier(diagram, group2,
            {{-1.0, -0.75, 0.0}, {0.85, 0.75, 0.0}, {-2.35, 0.5, 0.0},
                {1.0, -0.5, 0.0}});
        add_polyline(diagram, group2,
            {{-1.0, -0.75, 0.0}, {0.85, 0.75, 0.0}, {-2.35, 0.5, 0.0},
                {1.0, -0.5, 0.0}},
            {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!r", "$\\mathbf{p}_3$!!l",
                "$\\mathbf{p}_4$!!r"},
            stroke1);
      }
    }
    if (contains(tag, {"joins"})) {
      {
        auto group1 = add_group(diagram, {-3.0, 0.0, 0.0});
        add_qbezier(diagram, group1,
            {{-1.0, -0.75, 0.0}, {0.0, 0.75, 0.0}, {1.0, 0.0, 0.0}});
        add_polyline(diagram, group1,
            {{-1.0, -0.75, 0.0}, {0.0, 0.75, 0.0}, {1.0, 0.0, 0.0}},
            {"$\\mathbf{p}_{1,1}$!!r", "$\\mathbf{p}_{1,2}$!!l",
                "$\\mathbf{p}_{1,3}$!!l"},
            stroke1);
        add_point(diagram, group1, {0.0, 0.1875, 0.0}, "$\\mathbf{p}(u)$!!b");
      }
      {
        auto group2 = add_group(diagram, {1.0, 0.0, 0.0});
        add_bezier(diagram, group2,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}});
        add_point(diagram, group2, {-0.0375, 0.43, 0.0}, "$\\mathbf{p}(u)$!!b");
        add_polyline(diagram, group2,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}},
            {"$\\mathbf{p}_{1,1}$!!r", "$\\mathbf{p}_{1,2}$!!l",
                "$\\mathbf{p}_{1,3}$!!r", "$\\mathbf{p}_{1,4}$!!l"},
            stroke1);
      }
      {
        auto group3 = add_group(diagram, {-1, 0.0, 0.0});
        add_qbezier(diagram, group3,
            {{-1.0, 0.0, 0.0}, {0.0, -0.75, 0.0}, {1.0, 0.75, 0.0}});
        add_polyline(diagram, group3,
            {{-1.0, 0.0, 0.0}, {0.0, -0.75, 0.0}, {1.0, 0.75, 0.0}},
            {"$\\mathbf{p}_{2,1}$!!r", "$\\mathbf{p}_{2,2}$!!l",
                "$\\mathbf{p}_{2,3}$!!l"},
            stroke2);
      }
      {
        auto group4 = add_group(diagram, {3.0, 0.0, 0.0});
        add_bezier(diagram, group4,
            {{-1.0, 0.0, 0.0}, {-0.25, -0.75, 0.0}, {0.15, -0.75, 0.0},
                {1.0, 0.75, 0.0}});
        add_polyline(diagram, group4,
            {{-1.0, 0.0, 0.0}, {-0.25, -0.75, 0.0}, {0.15, -0.75, 0.0},
                {1.0, 0.75, 0.0}},
            {"$\\mathbf{p}_{2,1}$!!r", "$\\mathbf{p}_{2,2}$!!l",
                "$\\mathbf{p}_{2,3}$!!r", "$\\mathbf{p}_{2,4}$!!l"},
            stroke2);
      }
    }
    if (contains(tag, {"transforms"})) {
      {
        auto group1 = add_group(diagram, {-1.5, 0.0, 0.0});
        add_bezier(diagram, group1,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}});
        add_point(diagram, group1, {-0.0375, 0.43, 0.0}, "$\\mathbf{p}(u)$!!b");
        add_polyline(diagram, group1,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}},
            {"$\\mathbf{p}_1$!!l", "$\\mathbf{p}_2$!!l", "$\\mathbf{p}_3$!!r",
                "$\\mathbf{p}_4$!!r"},
            stroke1);
      }
      {
        auto group2 = add_group(
            diagram, {1.5, 0.0, 0.0}, drotation({0.0, 0.0, 305}));
        add_bezier(diagram, group2,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}});
        add_point(diagram, group2, {-0.0375, 0.43, 0.0},
            "$\\mathbf{p}(u,\\mathbf{F} \\cdot \\mathbf{q}_i) \\equiv$ \\\\ "
            "$\\quad \\equiv \\mathbf{F} \\cdot \\mathbf{p}(u,\\mathbf{q}_i)$!!bl");
        add_polyline(diagram, group2,
            {{-1.0, -0.75, 0.0}, {-0.35, 0.65, 0.0}, {0.25, 0.75, 0.0},
                {1.0, 0.0, 0.0}},
            {"$\\mathbf{F} \\cdot \\mathbf{p}_1$!!t",
                "$\\mathbf{F} \\cdot \\mathbf{p}_2$!!t",
                "$\\mathbf{F} \\cdot \\mathbf{p}_3$!!thr",
                "$\\mathbf{F} \\cdot \\mathbf{p}_4$!!r"},
            stroke1);
      }
    }
  }
  return diagram;
}

diagram_data subcurve_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto positions = vector<vec3f>{{-1.0, -1.0, 0.0}, {+1.0, -1.0, 0.0},
      {+1.0, +1.0, 0.0}, {-1.0, +1.0, 0.0}};
  auto lines     = vector<vec2i>{{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  auto [llines, lpositions]   = subdivide_lines(lines, positions, 1);
  auto [slines1, spositions1] = subdivide_bspline(lines, positions, 1);
  auto [slines2, spositions2] = subdivide_bspline(lines, positions, 2);
  auto [slines3, spositions3] = subdivide_bspline(lines, positions, 3);

  if (contains(tag, {"step"})) {
    {
      auto group1 = add_group(diagram, {-2.5, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "polyline!!t");
      add_lines(diagram, group1, lines, positions);
      add_points(diagram, group1, positions, {}, stroke1);
    }
    {
      auto group2 = add_group(diagram, {0, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "split!!t");
      add_lines(diagram, group2, llines, lpositions);
      add_points(diagram, group2, lpositions, {}, stroke1);
    }
    {
      auto group3 = add_group(diagram, {2.5, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "smooth!!t");
      add_lines(diagram, group3, slines1, spositions1);
      add_points(diagram, group3, spositions1, {}, stroke1);
    }
  }
  if (contains(tag, {"steps"})) {
    {
      auto group1 = add_group(diagram, {-2.5, 0.0, 0.0});
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "level 1!!t");
      add_lines(diagram, group1, slines1, spositions1);
      add_lines(diagram, group1, lines, positions, {}, stroke3);
    }
    {
      auto group2 = add_group(diagram, {0, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 1.0, 0.0}, "level 2!!t");
      add_lines(diagram, group2, slines2, spositions2);
      add_lines(diagram, group2, lines, positions, {}, stroke3);
    }
    {
      auto group3 = add_group(diagram, {2.5, 0.0, 0.0});
      add_label(diagram, group3, {0.0, 1.0, 0.0}, "level 3!!t");
      add_lines(diagram, group3, slines3, spositions3);
      add_lines(diagram, group3, lines, positions, {}, stroke3);
    }
  }
  return diagram;
}

diagram_data subdiv_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto positions              = dconstants::cube_positions;
  auto quads                  = dconstants::cube_quads;
  auto [lquads, lpositions]   = subdivide_quads(quads, positions, 1);
  auto [squads1, spositions1] = subdivide_catmullclark_(quads, positions, 1);
  auto [squads2, spositions2] = subdivide_catmullclark_(quads, positions, 2);
  auto [squads3, spositions3] = subdivide_catmullclark_(quads, positions, 3);
  auto [uquads3, upositions3] = subdivide_catmullclark_(
      quads, positions, 3, true);

  if (contains(tag, {"correction"})) {
    {
      auto group1 = add_group(
          diagram, {-1.5, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.25, 0.0}, "corrected!!t");
      add_quads(diagram, group1, squads3, spositions3);
      add_quads(diagram, group1, quads, positions, {}, transparent, stroke3);
    }
    {
      auto group2 = add_group(
          diagram, {1.5, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group2, {0.0, 1.25, 0.0}, "uncorrected!!t");
      add_quads(diagram, group2, uquads3, upositions3);
      add_quads(diagram, group2, quads, positions, {}, transparent, stroke3);
    }
  }
  if (contains(tag, {"step"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.25, 0.0}, "mesh!!t");
      add_quads(diagram, group1, quads, positions);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group2, {0.0, 1.25, 0.0}, "split!!t");
      add_quads(diagram, group2, lquads, lpositions);
    }
    {
      auto group3 = add_group(
          diagram, {3.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group3, {0.0, 1.25, 0.0}, "smooth!!t");
      add_quads(diagram, group3, squads1, spositions1);
    }
  }
  if (contains(tag, {"steps"})) {
    {
      auto group1 = add_group(
          diagram, {-3.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group1, {0.0, 1.25, 0.0}, "level 1!!t");
      add_quads(diagram, group1, squads1, spositions1);
      add_quads(diagram, group1, quads, positions, {}, transparent, stroke3);
    }
    {
      auto group2 = add_group(
          diagram, {0.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group2, {0.0, 1.25, 0.0}, "level 2!!t");
      add_quads(diagram, group2, squads2, spositions2);
      add_quads(diagram, group2, quads, positions, {}, transparent, stroke3);
    }
    {
      auto group3 = add_group(
          diagram, {3.0, 0.0, 0.0}, drotation({15.0, 300.0, 0.0}));
      add_label(diagram, group3, {0.0, 1.25, 0.0}, "level 3!!t");
      add_quads(diagram, group3, squads3, spositions3);
      add_quads(diagram, group3, quads, positions, {}, transparent, stroke3);
    }
  }
  return diagram;
}

diagram_data integration_diagrams(const string& tag) {
  auto func  = [](float x) -> float { return 1.75f - 1.5f * x * x / 9; };
  auto curve = vector<vec3f>{};
  for (auto idx = 0; idx <= 32; idx++) {
    auto x = 3 * idx / (float)32;
    curve.push_back({x, func(x), 0});
  }
  auto poly = curve;
  poly.push_back({3, 0, 0});
  poly.push_back({0, 0, 0});

  auto quadrature = tag == "quadrature";
  auto step       = 3.0f / 4.0f;
  auto points     = vector<vec3f>{};
  auto ticks      = quadrature ? vector<float>{0.5f, 1.5f, 2.5f, 3.5f}
                               : vector<float>{2.75f, 1.75f, 0.6f, 3.4f};
  for (auto idx = 0; idx < 4; idx++) {
    auto x = ticks[idx] * step;
    points.push_back({x, func(x), 0});
  }

  auto scale = 0.5f;

  auto diagram = make_diagram();

  {
    auto group1 = add_group(diagram, {-3.85, -1.0, 0.0});
    add_polygon(diagram, group1, poly, {}, fill1, transparent);
    add_label(diagram, group1, {1.5, 0.75, 0}, "$\\int f(x)dx$");
    add_polyline(diagram, group1, curve);
    add_arrows(
        diagram, group1, {{0, 0, 0}, {0, 2, 0}, {0, 0, 0}, {3.25, 0, 0}});
    add_labels(diagram, group1, {{0, 0, 0}, {3, 0, 0}}, {"$a$!!b", "$b$!!b"});
  }
  {
    auto group2 = add_group(diagram, {0.4, -1.0, 0.0});
    add_polyline(diagram, group2, curve);
    add_arrows(
        diagram, group2, {{0, 0, 0}, {0, 2, 0}, {0, 0, 0}, {3.25, 0, 0}});
    add_labels(diagram, group2, {{0, 0, 0}, {3, 0, 0}}, {"$a$!!b", "$b$!!b"});
    add_label(diagram, group2, {0.0, 2.0, 0.0},
        quadrature ? "$\\int f(x)dx \\approx \\sum_i f(x_i) (b-a)/n$!!t"
                   : "$\\int f(x)dx \\approx (1/n) \\sum_i f(x_i) (b-a)$!!t");
    {
      for (auto idx = 0; idx < points.size(); idx++) {
        auto [x, y, z] = points[idx];
        add_line(diagram, group2, {x, -0.05f, 0}, {x, +0.05f, 0});
        add_label(diagram, group2, {x, 0, 0},
            "$x_" + std::to_string(idx + 1) + "$!!b");
        if (quadrature) {
          add_label(diagram, group2, {x, y, 0},
              "$f(x_" + std::to_string(idx + 1) + ")$!!tr");
        } else {
          add_label(diagram, group2, {3, y, 0},
              "$f(x_" + std::to_string(idx + 1) + ")$!!r");
        }
      }
    }
    {
      for (auto idx = 0; idx < points.size(); idx++) {
        auto [x, y, z] = points[idx];
        if (quadrature) {
          add_quadv(diagram, group2,
              {{x - step / 2, y, 0}, {x + step / 2, y, 0}, {x + step / 2, 0, 0},
                  {x - step / 2, 0, 0}});
        } else {
          add_quadv(diagram, group2,
              {{0, y, -idx * 0.015f}, {3, y, -idx * 0.015f},
                  {3, 0, -idx * 0.015f}, {0, 0, -idx * 0.015f}});
        }
      }
    }
  }
  return diagram;
}

static float eval_pi(int samples, int trial) {
  auto sum = 0.0f;
  auto rng = make_rng(87687161, (size_t)trial * 2 + 1);
  for (auto sample = 0; sample < samples; sample++) {
    auto x = rand1f(rng) * 2 - 1, y = rand1f(rng) * 2 - 1;
    auto fxy = (x * x + y * y < 1) ? 1 : 0;
    sum += fxy;
  }
  return 4 * sum / samples;
};
static float eval_I(int samples, int trial) {
  auto sum = 0.0f;
  auto rng = make_rng(87687161, (size_t)trial * 2 + 1);
  for (auto sample = 0; sample < samples; sample++) {
    auto x  = rand1f(rng);
    auto fx = 1 - 0.75f * x * x;
    sum += fx;
  }
  return sum / samples;
};
const auto ifunc = 0.75f;

diagram_data mcplot_diagrams(const string& tag) {
  auto diagram = make_diagram();

  if (contains(tag, {"mcfuncerror"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"$I$ error"},
             .xbounds = {100.0, 10000.0},
             .xlabels = {"100", "10000", "samples"},
             .xticks  = {100.0, 10000.0, 5000.0},
             .ybounds = {0.65, 0.85},
             .ylabels = {"$I$", "$I-0.1$", "$I+0.1$"},
             .yticks  = {0.75, 0.65, 0.85}});
    add_lineplot(
        diagram, plot, [](float x) { return eval_I((int)x, (int)x); },
        {100.0, 10000.0}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return ifunc + ifunc / sqrt((float)x); },
        {100.0, 10000.0}, {}, stroke3);
    add_lineplot(
        diagram, plot, [](float x) { return ifunc - ifunc / sqrt((float)x); },
        {100.0, 10000.0}, {}, stroke3);
  }
  if (contains(tag, {"mcfunceval"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"$I$ estimate"},
             .xbounds = {1.0, 100.0},
             .xlabels = {"1", "100", "trial"},
             .xticks  = {1.0, 100.0, 50.0},
             .ybounds = {0.65, 0.85},
             .ylabels = {"$I$", "$I-0.1$", "$I+0.1$"},
             .yticks  = {0.75, 0.65, 0.85}});
    add_lineplot(
        diagram, plot, [](float x) { return eval_I(100, (int)x); },
        {1.0, 100.0}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return eval_I(1000, (int)x); },
        {1.0, 100.0}, {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return eval_I(10000, (int)x); },
        {1.0, 100.0}, {}, stroke3);
  }
  if (contains(tag, {"mcpierror"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"$\\pi$ error"},
             .xbounds = {1000.0, 100000.0},
             .xlabels = {"1000", "100000", "samples"},
             .xticks  = {1000.0, 100000.0, 50000.0},
             .ybounds = {pif - 0.2f, pif + 0.2f},
             .ylabels = {"$\\pi$", "$\\pi-0.2$", "$\\pi+0.2$"},
             .yticks  = {pif, pif - 0.2f, pif + 0.2f}});
    add_lineplot(
        diagram, plot, [](float x) { return eval_pi((int)x, (int)x); },
        {1000.0, 100000.0}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return pif + pif / sqrt((float)x); },
        {1000.0, 100000.0}, {}, stroke3);
    add_lineplot(
        diagram, plot, [](float x) { return pif - pif / sqrt((float)x); },
        {1000.0, 100000.0}, {}, stroke3);
  }
  if (contains(tag, {"mcpieval"})) {
    auto group = add_group(diagram);
    auto plot  = add_plot(diagram, group,
         {.titles     = {"$\\pi$ estimate"},
             .xbounds = {1.0, 100.0},
             .xlabels = {"1", "100", "trial"},
             .xticks  = {1.0, 100.0, 50.0},
             .ybounds = {pif - 0.2f, pif + 0.2f},
             .ylabels = {"$\\pi$", "$\\pi-0.2$", "$\\pi+0.2$"},
             .yticks  = {pif, pif - 0.2f, pif + 0.2f}});
    add_lineplot(
        diagram, plot, [](float x) { return eval_pi(1000, (int)x); },
        {1.0, 100.0}, {}, stroke1);
    add_lineplot(
        diagram, plot, [](float x) { return eval_pi(10000, (int)x); },
        {1.0, 100.0}, {}, stroke2);
    add_lineplot(
        diagram, plot, [](float x) { return eval_pi(100000, (int)x); },
        {1.0, 100.0}, {}, stroke3);
  }
  return diagram;
}

diagram_data pisamples_diagrams(const string& tag) {
  auto samples = 32 * 32;
  auto points  = vector<vec2f>{};
  auto rng     = make_rng(87687161);
  auto sum     = 0.0f;
  for (auto _ : range(samples)) {
    auto point = rand2f(rng) * 2 - 1;
    sum += dot(point, point) <= 1 ? 1 : 0;
    points.push_back(point);
  }
  auto estimate = 4 * sum / (float)samples;
  auto error    = abs(estimate - pif);

  auto steps   = vec2i{32, 32};
  auto spoints = vector<vec2f>{};
  auto ssum    = 0.0f;
  for (auto i : range(steps.x)) {
    for (auto j : range(steps.y)) {
      auto point =
          vec2f{(i + rand1f(rng)) / steps.x, (j + rand1f(rng)) / steps.y} * 2 -
          1;
      ssum += dot(point, point) <= 1 ? 1 : 0;
      spoints.push_back(point);
    }
  }
  auto sestimate = 4 * ssum / (float)samples;
  auto serror    = abs(sestimate - pif);

  auto diagram = make_diagram();

  if (contains(tag, {"stratified"})) {
    {
      auto group1 = add_group(diagram, {-1.2, 0.0, 0.0});
      add_points(diagram, group1, points, {}, black, dthickness::thin);
      add_quad(diagram, group1, {0.0, 0.0, -0.01}, 1);
      add_disk(diagram, group1, {0, 0, 0}, 1, {}, fill3);
    }
    {
      auto group2 = add_group(diagram, {-4.0, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 0.3, 0.0},
          "samples: " + std::to_string(samples) + "!!r");
      add_label(diagram, group2, {0.0, 0.0, 0.0},
          "pi: " + std::to_string(estimate) + "!!r");
      add_label(diagram, group2, {0.0, -0.3, 0.0},
          "error: " + std::to_string(error) + "!!r");
    }
    {
      auto group3 = add_group(diagram, {1.1, 0.0, 0.0});
      add_points(diagram, group3, spoints, {}, black, dthickness::thin);
      add_quad(diagram, group3, {0.0, 0.0, -0.01}, 1);
      add_disk(diagram, group3, {0, 0, 0}, 1, {}, fill3);
    }
    {
      auto group4 = add_group(diagram, {2.4, 0.0, 0.0});
      add_label(diagram, group4, {0.0, 0.3, 0.0},
          "samples: " + std::to_string(samples) + "!!r");
      add_label(diagram, group4, {0.0, 0.0, 0.0},
          "pi: " + std::to_string(sestimate) + "!!r");
      add_label(diagram, group4, {0.0, -0.3, 0.0},
          "error: " + std::to_string(serror) + "!!r");
    }
  }
  if (contains(tag, {"uniform"})) {
    {
      auto group1 = add_group(diagram, {-1.1, 0.0, 0.0});
      add_points(diagram, group1, points, {}, black, dthickness::thin);
      add_quad(diagram, group1, {0.0, 0.0, -0.01}, 1);
      add_disk(diagram, group1, {0, 0, 0}, 1, {}, fill3);
    }
    {
      auto group2 = add_group(diagram, {0.1, 0.0, 0.0});
      add_label(diagram, group2, {0.0, 0.3, 0.0},
          "samples: " + std::to_string(samples) + "!!r");
      add_label(diagram, group2, {0.0, 0.0, 0.0},
          "pi: " + std::to_string(estimate) + "!!r");
      add_label(diagram, group2, {0.0, -0.3, 0.0},
          "error: " + std::to_string(error) + "!!r");
    }
  }
  return diagram;
}

diagram_data sampling_diagrams(const string& tag) {
  auto diagram = make_diagram();

  auto nsamples = 16 * 16;

  {
    {
      auto group1 = add_group(diagram, {-2.3, 0.0, 0.0});
      add_rpoints(
          diagram, group1, {0, 0, 0}, 1, nsamples, false, {}, stroke1, thin);
      add_label(diagram, group1, {0.0, 1.0, 0.0}, "uniform square!!t");
      add_quad(diagram, group1, {0, 0, 0}, 1, {}, transparent);
      add_grid(diagram, group1, {0, 0, 0}, 1, {4, 4}, {}, gray, thin);
    }
    if (contains(tag, {"disk"})) {
      {
        auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
        add_rdpoints(
            diagram, group2, {0, 0, 0}, 1, nsamples, false, {}, stroke2, thin);
        add_label(diagram, group2, {0.0, 1.0, 0.0}, "uniform disk!!t");
        add_disk(diagram, group2, {0, 0, 0}, 1, {}, transparent);
        add_cgridnu(diagram, group2, {0, 0, 0}, 1, {4, 4}, {}, gray, thin);
      }
      {
        auto group3 = add_group(diagram, {2.3, 0.0, 0.0});
        add_rdpointsnu(
            diagram, group3, {0, 0, 0}, 1, nsamples, false, {}, stroke2, thin);
        add_label(diagram, group3, {0.0, 1.0, 0.0}, "non-uniform disk!!t");
        add_disk(diagram, group3, {0, 0, 0}, 1, {}, transparent);
        add_cgrid(diagram, group3, {0, 0, 0}, 1, {4, 4}, {}, gray, thin);
      }
    }
    if (contains(tag, {"stratified"})) {
      {
        auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
        add_rpoints(
            diagram, group2, {0, 0, 0}, 1, nsamples, true, {}, stroke2, thin);
        add_label(diagram, group2, {0.0, 1.0, 0.0}, "stratified samples!!t");
        add_quad(diagram, group2, {0, 0, 0}, 1, {}, transparent);
        add_grid(diagram, group2, {0, 0, 0}, 1, {4, 4}, {}, gray, thin);
      }
      {
        auto group3 = add_group(diagram, {2.3, 0.0, 0.0});
        add_rpoints(
            diagram, group3, {0, 0, 0}, 1, nsamples, true, {}, stroke2, thin);
        add_label(diagram, group3, {0.0, 1.0, 0.0}, "stratified grid!!t");
        add_quad(diagram, group3, {0, 0, 0}, 1, {}, transparent);
        add_grid(diagram, group3, {0, 0, 0}, 1, {16, 16}, {}, gray, thin);
      }
    }
    if (contains(tag, {"triangle"})) {
      {
        auto group2 = add_group(diagram, {0.0, 0.0, 0.0});
        add_rtpoints(
            diagram, group2, {0, 0, 0}, 1, nsamples, false, {}, stroke2, thin);
        add_label(diagram, group2, {0.0, 1.0, 0.0}, "triangle samples!!t");
        add_triangle(diagram, group2, {0, 0, 0}, 1, {}, transparent);
        add_tgrid(diagram, group2, {0, 0, 0}, 1, {4, 4}, {}, gray, thin);
      }
      {
        auto vertices = vector<vec3f>{
            {-1.0, -1.0, 0.0}, {+1.0, -1.0, 0.0}, {+1.0, +1.0, 0.0}};
        auto group3 = add_group(diagram, {2.3, 0.0, 0.0});
        add_rtpoints(
            diagram, group3, vertices, nsamples, false, {}, stroke3, thin);
        add_label(diagram, group3, {0.0, 1.0, 0.0}, "triangle samples!!t");
        add_trianglev(diagram, group3, vertices, {}, transparent);
        add_tgridv(diagram, group3, vertices, {4, 4}, {}, gray, thin);
      }
    }
  }
  return diagram;
}

}  // namespace yocto
