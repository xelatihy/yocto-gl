#include <yocto/yocto_cli.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>

using namespace yocto;

#ifdef _WIN32
#undef near
#undef far
#undef max
#endif

#include <yocto/yocto_color.h>

#include "yshade_sculpt.h"
#include "yshade_sculpt_algorithms.h"

// Application state
struct shade_sculpt_state {
  // loading parameters
  string shapename   = "";
  string filename    = "shape.ply";
  string imagename   = "";
  string outname     = "out.ply";
  string texturename = "";
  string name        = "";

  // scene
  sceneio_scene ioscene  = sceneio_scene{};
  camera_handle iocamera = invalid_handle;

  // options
  shade_params drawgl_prms = {};

  // rendering state
  shade_scene glscene = {};
};

struct sculpt_params {
  // brush type
  brush_type type = brush_type::gaussian;

  // intersection
  ray3f camera_ray = {};
  // sceneio_instance * shape_instance   = nullptr;
  shape_bvh          bvh          = {};
  shape_intersection intersection = {};

  // hash grid
  hash_grid grid = {};

  // stroke
  vec3f locked_position = {};
  vec2f locked_uv       = {};
  bool  lock            = false;
  bool  continuous      = false;
  bool  symmetric       = false;
  axes  symmetric_axis  = axes::x;

  // brush
  float               radius      = 0.35f;
  float               strength    = 1.0f;
  bool                negative    = false;
  vector<float>       opacity     = {};
  bool                saturation  = false;
  vector<vector<int>> adjacencies = {};

  // stroke parameterization
  vector<int>     stroke_sampling           = {};
  vector<int>     symmetric_stroke_sampling = {};
  geodesic_solver solver                    = {};
  vector<vec2f>   coords                    = {};
  image_data      tex_image                 = {};
  vector<vec3f>   old_positions             = {};
  vector<vec3f>   old_normals               = {};
};

// Initialize all sculpting parameters.
void init_sculpt_tool(
    sculpt_params &params, shape_data &shape, const scene_texture &texture) {
  // save positions
  params.old_positions = shape.positions;

  // create bvh structure
  params.bvh = make_triangles_bvh(
      shape.triangles, shape.positions, shape.radius);

  // create an hash grid
  params.grid = make_hash_grid(shape.positions, 0.05f);

  // init saturation buffer
  params.opacity = vector<float>(shape.positions.size(), 0.0f);

  // create geodesic distance graph (ONLY TRIANGLES MESHES!)
  auto adjacencies = face_adjacencies(shape.triangles);
  params.solver    = make_geodesic_solver(
      shape.triangles, adjacencies, shape.positions);
  params.adjacencies = vertex_adjacencies(shape.triangles, adjacencies);

  // init texture
  params.tex_image = texture;
}

shape_data make_circle(
    const vec3f &center, const mat3f &basis, float radius, int steps) {
  // 4 initial vertices
  auto  lines    = make_lines({1, 4});
  vec3f next_dir = basis.x;
  for (auto line : lines.lines) {
    auto v1             = line.x;
    auto v2             = line.y;
    lines.positions[v1] = center + next_dir * radius;
    lines.normals[v1]   = basis.z;
    next_dir            = cross(basis.z, next_dir);
    lines.positions[v2] = center + next_dir * radius;
    lines.normals[v2]   = basis.z;
  }

  // create polylines and fix lenght to radius
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto lins      = vector<vec2i>{};
  steps          = int(steps / 4);
  for (auto line : lines.lines) {
    auto v1  = line.x;
    auto v2  = line.y;
    auto pos = lines.positions[v1];
    positions.push_back(pos);
    normals.push_back(lines.normals[v1]);
    auto dir    = lines.positions[v2] - lines.positions[v1];
    auto lenght = length(dir) / steps;
    dir         = normalize(dir);
    for (int i = 0; i < steps; i++) {
      auto new_pos = pos + dir * lenght * (i + 1);
      auto new_dir = normalize(new_pos - center);
      new_pos      = center + new_dir * radius;
      positions.push_back(new_pos);
      normals.push_back(basis.z);
      lins.push_back({int(positions.size() - 2), int(positions.size() - 1)});
    }
  }

  // apply
  lines.positions = positions;
  lines.normals   = normals;
  lines.lines     = lins;
  return lines;
}

// To visualize mouse intersection on mesh
void view_pointer(shape_data &shape, shade_shape &glshape,
    shape_intersection intersection, float radius, int definition,
    brush_type &type) {
  if (intersection.hit) {
    if (type == brush_type::gaussian) radius *= 0.5f;
    auto pos    = eval_position(shape, intersection.element, intersection.uv);
    auto nor    = eval_normal(shape, intersection.element, intersection.uv);
    auto basis  = basis_fromz(nor);
    auto circle = make_circle(pos, basis, radius, definition);
    auto positions        = vector<vec3f>();
    positions             = circle.positions;
    auto lines            = circle.lines;
    auto nor_line         = make_lines({1, 1});
    nor_line.positions[0] = pos;
    nor_line.positions[1] = pos + nor * 0.05f;
    positions.insert(
        positions.end(), nor_line.positions.begin(), nor_line.positions.end());
    lines.push_back({int(positions.size() - 2), int(positions.size() - 1)});
    set_positions(glshape, positions);
    set_lines(glshape, lines);
  }
}

// To make the stroke sampling (position, normal) following the mouse
vector<pair<vec3f, vec3f>> stroke(sculpt_params &params, scene_shape &shape,
    const vec2f &mouse_uv, shade_camera &glcamera) {
  // eval current intersection
  auto  pairs = vector<pair<vec3f, vec3f>>{};
  auto &inter = params.intersection;
  auto  pos   = eval_position(
      shape.triangles, shape.positions, {inter.element, inter.uv});
  auto nor = eval_normal(
      shape.triangles, shape.normals, {inter.element, inter.uv});
  float delta_pos   = distance(pos, params.locked_position);
  float stroke_dist = params.radius * 0.2f;

  // handle continuous stroke
  if (params.continuous && (delta_pos < stroke_dist)) {
    pair<vec3f, vec3f> pair;
    pair.first  = pos;
    pair.second = nor;
    pairs.push_back(pair);
    params.stroke_sampling.push_back(
        closest_vertex(shape.triangles, inter.uv, inter.element));
    return pairs;
  }

  // handle first stroke intersection
  if (!params.lock) {
    params.locked_position = pos;
    params.locked_uv       = mouse_uv;
    params.lock            = true;
    return pairs;
  }

  float delta_uv = distance(mouse_uv, params.locked_uv);
  int   steps    = int(delta_pos / stroke_dist);
  if (steps == 0) return pairs;
  pairs           = vector<pair<vec3f, vec3f>>(steps);
  float stroke_uv = delta_uv * stroke_dist / delta_pos;
  auto  mouse_dir = normalize(mouse_uv - params.locked_uv);
  for (int step = 0; step < steps; step++) {
    params.locked_uv += stroke_uv * mouse_dir;
    auto ray = camera_ray(glcamera.frame, glcamera.lens, glcamera.aspect,
        glcamera.film, params.locked_uv);
    inter    = intersect_triangles_bvh(
        params.bvh, shape.triangles, shape.positions, ray, false);
    if (!inter.hit) continue;
    pos = eval_position(shape, inter.element, inter.uv);
    nor = eval_normal(shape, inter.element, inter.uv);
    params.stroke_sampling.push_back(
        closest_vertex(shape.triangles, inter.uv, inter.element));
    auto pair              = std::pair<vec3f, vec3f>{pos, nor};
    params.locked_position = pos;
    pairs.push_back(pair);
  }
  return pairs;
}

// End stroke settings
void end_stroke(sculpt_params &params, scene_shape &shape) {
  params.lock = false;
  std::fill(params.opacity.begin(), params.opacity.end(), 0.0f);
  params.stroke_sampling.clear();
  params.coords.clear();
  params.symmetric_stroke_sampling.clear();
  params.old_positions = shape.positions;
  params.old_normals   = shape.normals;
}

// Compute gaussian function
float gaussian_distribution(vec3f origin, vec3f position, float standard_dev,
    float scale_factor, float strength, float radius) {
  float scaled_strength = strength / ((((radius - 0.1f) * 0.5f) / 0.7f) + 0.2f);
  float N               = 1.0f / (((standard_dev * scaled_strength) *
                        (standard_dev * scaled_strength) *
                        (standard_dev * scaled_strength)) *
                       sqrt((2.0f * pi) * (2.0f * pi) * (2.0f * pi)));
  float dx              = (origin.x - position.x) * scale_factor;
  float dy              = (origin.y - position.y) * scale_factor;
  float dz              = (origin.z - position.z) * scale_factor;
  float E               = ((dx * dx) + (dy * dy) + (dz * dz)) /
            (2.0f * standard_dev * standard_dev);
  return N * yocto::exp(-E);
}

// Change positions, normals, bounding volume hierarchy and hash grid
void apply_brush(shape_data &shape, vector<vec3f> &positions, shape_bvh &tree,
    hash_grid &grid) {
  shape.positions = positions;
  triangles_normals(shape.normals, shape.triangles, shape.positions);
  update_triangles_bvh(tree, shape.triangles, shape.positions);
  grid = make_hash_grid(positions, grid.cell_size);
}

// To apply brush on intersected points' neighbors
void brush(sculpt_params &params, scene_shape &shape,
    vector<pair<vec3f, vec3f>> &pairs) {
  if (pairs.empty()) return;
  auto &positions = shape.positions;
  auto  neighbors = vector<int>{};

  // for a correct gaussian distribution
  float scale_factor = 3.5f / params.radius;

  for (int p = 0; p < pairs.size(); p++) {
    find_neighbors(params.grid, neighbors, pairs[p].first, params.radius);
    auto center = pairs[p].first;
    auto normal = pairs[p].second;
    if (params.negative) normal = -normal;
    if (params.saturation) {
      float max_height = gaussian_distribution(
          center, center, 0.7f, scale_factor, params.strength, params.radius);
      for (int v = 0; v < neighbors.size(); v++) {
        float opacity_val = params.opacity[neighbors[v]];
        if (opacity_val == 1.0f) continue;
        float gauss_height = gaussian_distribution(center,
            positions[neighbors[v]], 0.7f, scale_factor, params.strength,
            params.radius);
        float gauss_ratio  = gauss_height / max_height;
        if ((opacity_val + gauss_ratio) > 1.0f) {
          gauss_height = max_height * (1.0f - opacity_val);
          gauss_ratio  = 1.0f - opacity_val;
        }
        params.opacity[neighbors[v]] += gauss_ratio;
        positions[neighbors[v]] += normal * gauss_height;
      }
    } else {
      for (int v = 0; v < neighbors.size(); v++) {
        float gauss_height = gaussian_distribution(center,
            positions[neighbors[v]], 0.7f, scale_factor, params.strength,
            params.radius);
        positions[neighbors[v]] += normal * gauss_height;
      }
    }
    neighbors.clear();
  }

  apply_brush(shape, positions, params.bvh, params.grid);
}

// Compute texture values through the parameterization
void texture_brush(vector<int> &vertices, image_data &texture,
    vector<vec2f> &coords, sculpt_params &params, scene_shape &shape,
    vector<vec3f> positions, vector<vec3f> normals) {
  if (vertices.empty()) return;
  if (texture.pixelsf.empty() && texture.pixelsb.empty()) return;

  auto scale_factor = 3.5f / params.radius;
  auto max_height   = gaussian_distribution(
      zero3f, zero3f, 0.7f, scale_factor, params.strength, params.radius);

  for (auto i : vertices) {
    auto uv     = coords[i];
    auto color  = eval_image(texture, uv);
    auto height = length(color);
    auto normal = normals[i];
    if (params.negative) normal = -normal;
    height *= max_height;
    positions[i] += normal * height;
  }

  apply_brush(shape, positions, params.bvh, params.grid);
}

// Cotangent operator
float cotan(vec3f &a, vec3f &b) { return dot(a, b) / length(cross(a, b)); }

// Compute edge cotangents weights
float laplacian_weight(vector<vec3f> &positions,
    vector<vector<int>> &adjacencies, int node, int neighbor) {
  auto num_neighbors = int(adjacencies[node].size());

  int ind = -1;
  for (int i = 0; i < adjacencies[node].size(); i++) {
    if (adjacencies[node][i] == neighbor) {
      ind = i;
      break;
    }
  }

  auto prev = adjacencies[node][(ind - 1) >= 0 ? ind - 1 : num_neighbors - 1];
  auto next = adjacencies[node][(ind + 1) % num_neighbors];

  auto v1 = positions[node] - positions[prev];
  auto v2 = positions[neighbor] - positions[prev];
  auto v3 = positions[node] - positions[next];
  auto v4 = positions[neighbor] - positions[next];

  float cot_alpha = cotan(v1, v2);
  float cot_beta  = cotan(v3, v4);
  float weight    = (cot_alpha + cot_beta) / 2;

  float cotan_max = yocto::cos(flt_min) / yocto::sin(flt_min);
  weight          = clamp(weight, cotan_max, -cotan_max);
  return weight;
}

// Smooth brush with Laplace Operator Discretization and Cotangents Weights
void smooth(geodesic_solver &solver, vector<int> &stroke_sampling,
    vector<vec3f> &positions, sculpt_params &params, scene_shape &shape) {
  if (stroke_sampling.empty()) return;

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : stroke_sampling) distances[sample] = 0.0f;

  int  current_node = -1;
  auto neighbors    = vector<int>{};
  auto weights      = vector<float>{};
  auto update       = [&](int node, int neighbor, float new_distance) {
    if (current_node == -1) current_node = node;
    if (node != current_node) {
      vec3f sum1 = zero3f;
      float sum2 = 0.0f;
      for (int i = 0; i < neighbors.size(); i++) {
        sum1 += positions[neighbors[i]] * weights[i];
        sum2 += weights[i];
      }
      positions[current_node] += 0.5f *
                                 (((sum1 / sum2) - positions[current_node]));
      current_node = node;
      neighbors.clear();
      weights.clear();
    }
    neighbors.push_back(neighbor);
    weights.push_back(
        laplacian_weight(positions, params.adjacencies, node, neighbor));
  };
  dijkstra(solver, stroke_sampling, distances, params.radius, update);

  apply_brush(shape, positions, params.bvh, params.grid);

  stroke_sampling.clear();
}

// TODO(fabio): move this function to math
static frame3f camera_frame(float lens, float aspect, float film = 0.036) {
  auto camera_dir  = normalize(vec3f{0, 0.5, 1});
  auto bbox_radius = 2.0f;
  auto camera_dist = bbox_radius * lens / (film / aspect);
  return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
}

static void convert_scene(scene_scene &scene, const scene_shape &ioshape_,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 5};
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);

  // init scene
  scene = {};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto &pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto &pos : ioshape.positions) pos -= center(bbox);
  for (auto &pos : ioshape.positions) pos /= max(size(bbox));

  // camera
  if (progress_cb) progress_cb("create camera", progress.x++, progress.y);
  auto &camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  if (progress_cb) progress_cb("create material", progress.x++, progress.y);
  auto &shape_material  = scene.materials.emplace_back();
  shape_material.type   = material_type::matte;
  shape_material.color  = {0.78f, 0.31f, 0.23f};
  auto &cursor_material = scene.materials.emplace_back();
  cursor_material.type  = material_type::matte;

  // shapes
  if (progress_cb) progress_cb("create shape", progress.x++, progress.y);
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(
      make_circle({0, 0, 0}, basis_fromz({0, 0, 1}), 1, 32));

  // instances
  if (progress_cb) progress_cb("create instance", progress.x++, progress.y);
  auto &shape_instance     = scene.instances.emplace_back();
  shape_instance.shape     = 0;
  shape_instance.material  = 0;
  auto &cursor_instance    = scene.instances.emplace_back();
  cursor_instance.shape    = 1;
  cursor_instance.material = 1;

  // done
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);
}

static void init_glscene(shade_scene &glscene, const sceneio_scene &ioscene,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene.cameras.size() + (int)ioscene.materials.size() +
             (int)ioscene.textures.size() + (int)ioscene.shapes.size() +
             (int)ioscene.instances.size()};

  // init scene
  init_scene(glscene);

  // camera
  for (auto &iocamera : ioscene.cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto &camera = glscene.cameras.at(add_camera(glscene));
    set_frame(camera, iocamera.frame);
    set_lens(camera, iocamera.lens, iocamera.aspect, iocamera.film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  for (auto &iotexture : ioscene.textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto  handle    = add_texture(glscene);
    auto &gltexture = glscene.textures[handle];
    if (!iotexture.pixelsf.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsf);
    } else if (!iotexture.pixelsb.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsb);
    }
  }

  // material
  for (auto &iomaterial : ioscene.materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto  handle     = add_material(glscene);
    auto &glmaterial = glscene.materials[handle];
    set_emission(glmaterial, iomaterial.emission, iomaterial.emission_tex);
    set_opacity(glmaterial, iomaterial.opacity, invalid_handle);
    set_normalmap(glmaterial, iomaterial.normal_tex);
    switch (iomaterial.type) {
      case material_type::matte: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, 0, invalid_handle);
      } break;
      case material_type::plastic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      case material_type::metal: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 1, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
      case material_type::metallic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, iomaterial.metallic, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      default: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
    }
  }

  // shapes
  for (auto &ioshape : ioscene.shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    add_shape(glscene, ioshape.points, ioshape.lines, ioshape.triangles,
        ioshape.quads, ioshape.positions, ioshape.normals, ioshape.texcoords,
        ioshape.colors);
  }

  // shapes
  for (auto &ioinstance : ioscene.instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto  handle     = add_instance(glscene);
    auto &glinstance = glscene.instances[handle];
    set_frame(glinstance, ioinstance.frame);
    set_shape(glinstance, ioinstance.shape);
    set_material(glinstance, ioinstance.material);
  }

  // environments
  for (auto &ioenvironment : ioscene.environments) {
    auto  handle        = add_environment(glscene);
    auto &glenvironment = glscene.environments[handle];
    set_frame(glenvironment, ioenvironment.frame);
    set_emission(
        glenvironment, ioenvironment.emission, ioenvironment.emission_tex);
  }

  // init environments
  init_environments(glscene);

  // done
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);
}

int run_shade_sculpt(const shade_sculpt_params &params_) {
  // initialize app
  auto app = shade_sculpt_state();

  // copy command line
  app.filename    = params_.shape;
  app.texturename = params_.texture;
  app.imagename   = params_.texture;
  app.shapename   = params_.shape;

  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  print_progress("load shape", 0, 1);
  if (!load_shape(app.filename, ioshape, ioerror)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads.clear();
  }
  print_progress("load shape", 1, 1);

  // loading texture
  auto iotexture = scene_texture{};
  if (!app.texturename.empty()) {
    print_progress("load texture", 0, 1);
    if (!load_texture(app.texturename, iotexture, ioerror))
      print_fatal(ioerror);
    print_progress("load texture", 1, 1);
  }

  // setup app
  convert_scene(app.ioscene, ioshape, print_progress);

  // sculpt params
  auto params = sculpt_params{};
  init_sculpt_tool(params, app.ioscene.shapes.front(), iotexture);

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&app](gui_window *win, const gui_input &input) {
    init_glscene(app.glscene, app.ioscene, print_progress);
  };
  callbacks.clear_cb = [&app](gui_window *win, const gui_input &input) {
    clear_scene(app.glscene);
  };
  callbacks.draw_cb = [&app](gui_window *win, const gui_input &input) {
    draw_scene(app.glscene, app.glscene.cameras.at(0),
        input.framebuffer_viewport, app.drawgl_prms);
  };
  callbacks.widgets_cb = [&app, &params](
                             gui_window *win, const gui_input &input) {
    auto &glparams = app.drawgl_prms;
    draw_slider(win, "resolution", glparams.resolution, 0, 4096);
    draw_checkbox(win, "wireframe", glparams.wireframe);
    draw_combobox(
        win, "lighting", (int &)glparams.lighting, shade_lighting_names);
    continue_line(win);
    draw_checkbox(win, "double sided", glparams.double_sided);
    draw_slider(win, "exposure", glparams.exposure, -10, 10);
    draw_slider(win, "gamma", glparams.gamma, 0.1f, 4);
    draw_slider(win, "near", glparams.near, 0.01f, 1.0f);
    draw_slider(win, "far", glparams.far, 1000.0f, 10000.0f);
    draw_label(win, "", "");
    draw_label(win, "", "sculpt params");
    draw_combobox(win, "brush type", (int &)params.type, brushes_names);
    if (params.type == brush_type::gaussian) {
      if (params.strength < 0.8f || params.strength > 1.5f)
        params.strength = 1.0f;
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", params.negative);
      draw_checkbox(win, "continuous", params.continuous);
      draw_checkbox(win, "saturation", params.saturation);
      draw_checkbox(win, "symmetric", params.symmetric);
      if (params.symmetric)
        draw_combobox(win, "symmetric axis", (int &)params.symmetric_axis,
            symmetric_axes);
    } else if (params.type == brush_type::texture) {
      if (params.strength < 0.8f || params.strength > 1.5f)
        params.strength = 1.0f;
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", params.negative);
      draw_checkbox(win, "symmetric", params.symmetric);
      if (params.symmetric)
        draw_combobox(win, "symmetric axis", (int &)params.symmetric_axis,
            symmetric_axes);
    } else if (params.type == brush_type::smooth) {
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 0.1f, 1.0f);
      draw_checkbox(win, "symmetric", params.symmetric);
      if (params.symmetric)
        draw_combobox(win, "symmetric axis", (int &)params.symmetric_axis,
            symmetric_axes);
    }
  };
  callbacks.update_cb = [](gui_window *win, const gui_input &input) {
    // pass
  };
  callbacks.uiupdate_cb = [&app, &params](
                              gui_window *win, const gui_input &input) {
    // intersect mouse position and shape
    auto mouse_uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
        input.mouse_pos.y / float(input.window_size.y)};

    auto &shape         = app.ioscene.shapes.at(0);
    auto &glcamera      = app.glscene.cameras.at(0);
    auto &glshape       = app.glscene.shapes.at(0);
    params.camera_ray   = camera_ray(glcamera.frame, glcamera.lens,
        glcamera.aspect, glcamera.film, mouse_uv);
    params.intersection = intersect_triangles_bvh(
        params.bvh, shape.triangles, shape.positions, params.camera_ray, false);
    if (params.intersection.hit) {
      app.glscene.instances.back().hidden = false;
      view_pointer(shape, app.glscene.shapes.back(), params.intersection,
          params.radius, 20, params.type);
    } else {
      app.glscene.instances.back().hidden = true;
    }

    auto isec = params.intersection;
    // sculpting
    if (input.mouse_left && isec.hit && !input.modifier_ctrl &&
        (isec.uv.x >= 0 && isec.uv.x < 1) &&
        (isec.uv.y >= 0 && isec.uv.y < 1)) {
      auto pairs = stroke(params, shape, mouse_uv, app.glscene.cameras.at(0));
      vector<int> vertices;
      if (params.type == brush_type::gaussian) {
        brush(params, shape, pairs);
      } else if (params.type == brush_type::smooth) {
        smooth(params.solver, params.stroke_sampling, shape.positions, params,
            shape);
      } else if (params.type == brush_type::texture && !pairs.empty()) {
        vertices = stroke_parameterization(params.solver, params.coords,
            params.stroke_sampling, params.old_positions, params.old_normals,
            params.radius);
        texture_brush(vertices, params.tex_image, params.coords, params, shape,
            params.old_positions, params.old_normals);
      }
      set_positions(glshape, shape.positions);
      set_normals(glshape, shape.normals);
      if (params.symmetric) {
        pairs = symmetric_stroke(pairs, shape, params.bvh,
            params.symmetric_stroke_sampling, params.symmetric_axis);
        if (params.type == brush_type::gaussian) {
          brush(params, shape, pairs);
        } else if (params.type == brush_type::smooth) {
          smooth(params.solver, params.symmetric_stroke_sampling,
              shape.positions, params, shape);
        } else if (params.type == brush_type::texture && !pairs.empty()) {
          vertices = stroke_parameterization(params.solver, params.coords,
              params.symmetric_stroke_sampling, params.old_positions,
              params.old_normals, params.radius);
          texture_brush(vertices, params.tex_image, params.coords, params,
              shape, shape.positions, shape.normals);
        }
        set_positions(glshape, shape.positions);
        set_normals(glshape, shape.normals);
      }
    } else {
      end_stroke(params, shape);
    }
    if (input.modifier_ctrl && !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      auto &glcamera                           = app.glscene.cameras.at(0);
      std::tie(glcamera.frame, glcamera.focus) = camera_turntable(
          glcamera.frame, glcamera.focus, rotate, dolly, -pan);
    }
  };

  // run ui
  run_ui({1900, 876}, "yshade", callbacks);

  // done
  return 0;
}
