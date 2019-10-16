#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
using namespace yocto;

struct app_state {
  // Callbacks available for user to build its own behaviors
  std::function<void(app_state&)>                         init;
  std::function<void(app_state&, int, bool)>              key_callback;
  std::function<void(app_state&, int, vec2f, int, float)> click_callback;
  std::function<void(app_state&, const opengl_window&)>   draw_glwidgets;

  // Geometry data
  yocto_shape shape;

  // OpenGL data
  opengl_scene        scene          = {};
  draw_glscene_params opengl_options = {};

  // Interaction data
  float        time       = 0;
  bool         show_edges = false;
  yocto_camera camera;
  float        camera_focus;
  bvh_tree     bvh;

  // Internal handles
  int glshape_id, glpoints_id, glvector_field_id, gledges_id, glpolyline_id;
  opengl_shape& glshape() { return scene.shapes[glshape_id]; }
  opengl_shape& glpoints() { return scene.shapes[glpoints_id]; }
  opengl_shape& glvector_field() { return scene.shapes[glvector_field_id]; }
  opengl_shape& gledges() { return scene.shapes[gledges_id]; }
  opengl_shape& glpolyline() { return scene.shapes[glpolyline_id]; }
};

// @Issue: Maybe this in yocto_opengl.h?
void delete_glshape(opengl_shape& glshape) {
  delete_glarraybuffer(glshape.positions);
  delete_glarraybuffer(glshape.normals);
  delete_glarraybuffer(glshape.texcoords);
  delete_glarraybuffer(glshape.colors);
  delete_glarraybuffer(glshape.tangentsps);
  delete_glelementbuffer(glshape.points);
  delete_glelementbuffer(glshape.lines);
  delete_glelementbuffer(glshape.triangles);
  delete_glelementbuffer(glshape.quads);
  delete_glelementbuffer(glshape.edges);
}

void update_glshape(app_state& app) {
  // @Issue: This app is specialized for a model that is a triangle mesh.
  //    Loading a generic shape is unsafe, maybe we should load only
  //    triangle meshes here...

  auto& glshape = app.glshape();
  auto& shape   = app.shape;
  delete_glshape(glshape);
  if (shape.quadspos.empty()) {
    if (!shape.positions.empty())
      init_glarraybuffer(glshape.positions, shape.positions, false);
    if (!shape.normals.empty())
      init_glarraybuffer(glshape.normals, shape.normals, false);
    if (!shape.texcoords.empty())
      init_glarraybuffer(glshape.texcoords, shape.texcoords, false);
    if (!shape.colors.empty())
      init_glarraybuffer(glshape.colors, shape.colors, false);
    if (!shape.tangents.empty())
      init_glarraybuffer(glshape.tangentsps, shape.tangents, false);
    if (!shape.points.empty())
      init_glelementbuffer(glshape.points, shape.points, false);
    if (!shape.lines.empty())
      init_glelementbuffer(glshape.lines, shape.lines, false);
    if (!shape.triangles.empty())
      init_glelementbuffer(glshape.triangles, shape.triangles, false);
    if (!shape.quads.empty()) {
      auto triangles = quads_to_triangles(shape.quads);
      init_glelementbuffer(glshape.quads, triangles, false);
    }
  } else {
    auto quads     = vector<vec4i>{};
    auto positions = vector<vec3f>{};
    auto normals   = vector<vec3f>{};
    auto texcoords = vector<vec2f>{};
    split_facevarying(quads, positions, normals, texcoords, shape.quadspos,
        shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
        shape.texcoords);
    if (!positions.empty())
      init_glarraybuffer(glshape.positions, positions, false);
    if (!normals.empty()) init_glarraybuffer(glshape.normals, normals, false);
    if (!texcoords.empty())
      init_glarraybuffer(glshape.texcoords, texcoords, false);
    if (!quads.empty()) {
      auto triangles = quads_to_triangles(quads);
      init_glelementbuffer(glshape.quads, triangles, false);
    }
  }
}

void update_glpolyline(app_state& app, const vector<vec3f>& vertices) {
  auto& glshape = app.glpolyline();
  delete_glshape(glshape);
  if (vertices.size()) {
    auto elements = vector<vec2i>(vertices.size() - 1);
    for (int i = 0; i < elements.size(); i++) elements[i] = {i, i + 1};
    init_glarraybuffer(glshape.positions, vertices, false);
    init_glelementbuffer(glshape.lines, elements, false);
  }
}

void update_glpoints(app_state& app, const vector<vec3f>& points) {
  auto& glshape = app.glpoints();
  delete_glshape(glshape);
  if (points.size()) {
    auto elements = vector<int>(points.size());
    for (int i = 0; i < elements.size(); i++) elements[i] = i;
    init_glarraybuffer(glshape.positions, points, false);
    init_glarraybuffer(
        glshape.normals, vector<vec3f>(points.size(), {0, 0, 1}), false);
    init_glelementbuffer(glshape.points, elements, false);
  }
}

void update_glvector_field(
    app_state& app, const vector<vec3f>& vector_field, float scale = 0.01) {
  auto perface   = vector_field.size() == app.shape.triangles.size();
  auto pervertex = vector_field.size() == app.shape.positions.size();

  if (!perface && !pervertex) {
    throw std::runtime_error("input vector field has wrong size\n");
  }

  auto& glshape = app.glvector_field();
  delete_glshape(glshape);
  auto size = perface ? app.shape.triangles.size() : app.shape.positions.size();
  auto positions = vector<vec3f>(size * 2);

  // Per-face vector field
  if (perface) {
    for (int i = 0; i < app.shape.triangles.size(); i++) {
      auto x      = app.shape.positions[app.shape.triangles[i].x];
      auto y      = app.shape.positions[app.shape.triangles[i].y];
      auto z      = app.shape.positions[app.shape.triangles[i].z];
      auto normal = triangle_normal(x, y, z);
      normal *= 0.0001;
      auto center          = (x + y + z) / 3;
      auto from            = center + normal;
      auto to              = from + (scale * vector_field[i]) + normal;
      positions[i * 2]     = from;
      positions[i * 2 + 1] = to;
    }
  } else {
    for (int i = 0; i < app.shape.positions.size(); i++) {
      auto from            = app.shape.positions[i];
      auto to              = from + scale * vector_field[i];
      positions[i * 2]     = from;
      positions[i * 2 + 1] = to;
    }
  }
  init_glarraybuffer(glshape.positions, positions, false);

  auto elements = vector<vec2i>(size);
  for (int i = 0; i < elements.size(); i++) {
    elements[i] = {2 * i, 2 * i + 1};
  }
  init_glelementbuffer(glshape.lines, elements, false);
}

void update_gledges(app_state& app) {
  auto& glshape = app.gledges();
  delete_glshape(glshape);
  auto positions = app.shape.positions;
  for (int i = 0; i < positions.size(); i++) {
    positions[i] += app.shape.normals[i] * 0.0001;
  }
  init_glarraybuffer(glshape.positions, positions, false);

  auto elements = vector<vec2i>();
  elements.reserve(app.shape.triangles.size() * 3);
  for (int i = 0; i < app.shape.triangles.size(); i++) {
    for (int k = 0; k < 3; k++) {
      auto a = app.shape.triangles[i][k];
      auto b = app.shape.triangles[i][(k + 1) % 3];
      if (a < b) {
        elements.push_back({a, b});
      }
    }
  }
  init_glelementbuffer(glshape.lines, elements, false);
}

void update_glcamera(opengl_camera& glcamera, const yocto_camera& camera) {
  glcamera.frame  = camera.frame;
  glcamera.yfov   = camera_yfov(camera);
  glcamera.asepct = camera_aspect(camera);
  glcamera.near   = 0.001f;
  glcamera.far    = 10000;
}

void init_camera(app_state& app, const vec3f& from = vec3f{0, 0.5, 1.5},
    const vec3f& to = {0, 0, 0}) {
  app.camera              = yocto_camera{};
  auto up                 = vec3f{0, 1, 0};
  app.camera.lens         = 0.02f;
  app.camera.orthographic = false;
  app.camera.aperture     = 0;
  app.camera.frame        = lookat_frame(from, to, up);
  app.camera.film         = {0.036f, 0.015f};
  app.camera.focus        = length(to - from);
  app.camera_focus        = app.camera.focus;
}

void init_bvh(app_state& app) {
  make_triangles_bvh(app.bvh, app.shape.triangles, app.shape.positions,
      app.shape.radius, false, false);
}

void hide_edges(app_state& app) {
  app.show_edges                            = false;
  app.scene.instances[app.gledges_id].shape = -1;
}
void show_edges(app_state& app) {
  app.show_edges                            = true;
  app.scene.instances[app.gledges_id].shape = app.gledges_id;
}

void init_opengl_scene(app_state& app) {
  make_glscene(app.scene);
  update_glcamera(app.scene.cameras.emplace_back(), app.camera);

  auto shape_material      = opengl_material{};
  shape_material.diffuse   = {1, 0.2, 0};
  shape_material.roughness = 0.3;
  app.scene.materials.push_back(shape_material);

  // @Issue: Right now we're missing APIs to color things easily.
  auto lines_material      = opengl_material{};
  lines_material.emission  = {1, 1, 1};
  lines_material.roughness = 0.0;
  app.scene.materials.push_back(lines_material);

  // The model.
  app.glshape_id = app.scene.shapes.size();
  app.scene.shapes.push_back({});
  update_glshape(app);

  // The points.
  app.glpoints_id = app.scene.shapes.size();
  app.scene.shapes.push_back({});

  // The vector field.
  app.glvector_field_id = app.scene.shapes.size();
  app.scene.shapes.push_back({});

  // The edges.
  app.gledges_id = app.scene.shapes.size();
  app.scene.shapes.push_back({});
  update_gledges(app);

  // The polyline.
  app.glpolyline_id = app.scene.shapes.size();
  app.scene.shapes.push_back({});

  // Add instances.
  app.scene.instances = vector<opengl_instance>(5);
  for (int i = 0; i < app.scene.instances.size(); ++i) {
    app.scene.instances[i].shape    = i;
    app.scene.instances[i].material = i ? 1 : 0;
  }

  // Hide edges.
  if (not app.show_edges) hide_edges(app);

  // Add lights.
  app.scene.lights.push_back({{5, 5, 5}, {30, 30, 30}, 0});
  app.scene.lights.push_back({{-5, 5, 5}, {30, 30, 30}, 0});
  app.scene.lights.push_back({{0, 5, -5}, {30, 30, 30}, 0});
}

void clear(app_state& app) {
  for (int i = 0; i < app.scene.shapes.size(); i++) {
    if (i == app.glshape_id) continue;
    if (i == app.gledges_id) continue;
    delete_glshape(app.scene.shapes[i]);
  }
  delete_glarraybuffer(app.glshape().colors);
  init_glarraybuffer(app.glshape().colors,
      vector<vec4f>(app.shape.positions.size(), {1, 1, 1, 1}));
}

void key_callback(const opengl_window& win, int key, bool pressing) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  app.key_callback(app, key, pressing);
}

void click_callback(const opengl_window& win, bool left_click, bool press) {
  auto& app   = *(app_state*)get_gluser_pointer(win);
  auto  mouse = get_glmouse_pos_normalized(win, false);

  // Ray trace camera ray.
  if (!left_click && press) {
    auto  ray = eval_camera(app.camera, mouse, {0.5, 0.5});
    int   face;
    vec2f uv;
    float distance;
    auto  hit = intersect_triangles_bvh(app.bvh, app.shape.triangles,
        app.shape.positions, ray, face, uv, distance);

    if (hit) {
      auto uvw = vec3f{uv.x, uv.y, 1 - uv.x - uv.y};
      int  k   = 0;
      if (uvw.x > uvw.y && uvw.x > uvw.z) k = 1;
      if (uvw.y > uvw.x && uvw.y > uvw.z) k = 2;
      auto vertex = app.shape.triangles[face][k];
      app.click_callback(app, face, uv, vertex, distance);
    }
  }
}

void scroll_callback(const opengl_window& win, float yoffset) {
  auto& app  = *(app_state*)get_gluser_pointer(win);
  float zoom = yoffset > 0 ? 0.1 : -0.1;
  update_turntable(app.camera.frame, app.camera.focus, zero2f, zoom, zero2f);
  update_glcamera(app.scene.cameras[0], app.camera);
}

void draw(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  draw_glscene(
      app.scene, get_glframebuffer_viewport(win, false), app.opengl_options);
  begin_glwidgets(win);
  app.draw_glwidgets(app, win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void run_app(app_state& app) {
  // Init window.
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimshproc", &app, draw);
  init_opengl_scene(app);

  set_click_glcallback(win, click_callback);
  set_scroll_glcallback(win, scroll_callback);
  set_key_glcallback(win, key_callback);

  init_glwidgets(win);

  auto mouse_pos = zero2f, last_pos = zero2f;

  // Render loop.
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto alt_down       = get_glalt_key(win);
    auto shift_down     = get_glshift_key(win);
    auto widgets_active = get_glwidgets_active(win);

    // Handle mouse and keyboard for navigation.
    if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
      auto& camera = app.camera;
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
      rotate.y = -rotate.y;
      pan.x    = -pan.x;
      update_turntable(camera.frame, app.camera.focus, rotate, dolly, pan);
      update_glcamera(app.scene.cameras[0], camera);
    }

    draw(win);

    // Event handling.
    process_glevents(win);
  }

  delete_glwindow(win);
}

void yimshproc(const string&                                  input_filename,
    std::function<void(app_state&)>                           init,
    std::function<void(app_state&, int, bool)>                key_callback,
    std::function<void(app_state&, int, vec2f, int, float)>   click_callback,
    std::function<void(app_state&, const opengl_window& win)> draw_glwidgets) {
  auto app = app_state{};

  // init shape
  load_shape(input_filename, app.shape.points, app.shape.lines,
      app.shape.triangles, app.shape.quads, app.shape.positions,
      app.shape.normals, app.shape.texcoords, app.shape.colors,
      app.shape.radius);
  init_bvh(app);
  init_camera(app);

  app.init           = init;
  app.key_callback   = key_callback;
  app.click_callback = click_callback;
  app.draw_glwidgets = draw_glwidgets;  // @Issue: win not needed for widgets

  app.init(app);
  run_app(app);
}
