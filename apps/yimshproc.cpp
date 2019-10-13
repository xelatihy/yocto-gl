#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <GLFW/glfw3.h>

struct app_state {
  string input_filename  = "model.obj";
  string output_filename = "model.obj";

  opengl_scene        scene          = {};
  draw_glscene_params opengl_options = {};

  float       time             = 0;
  vector<int> vertex_selection = {};

  yocto_shape         shape;
  vector<vec3i>       face_adjacency;
  vector<vector<int>> vertex_adjacency;
  geodesic_solver     solver;
  vector<float>       scalar_field;
  vector<vec3f>       vector_field;
  float               vector_field_scale = 1;
  bool                show_edges         = false;

  yocto_camera camera;
  float        camera_focus;
  bvh_tree     bvh;

  opengl_material& shape_material() { return scene.materials[0]; }
  opengl_material& lines_material() { return scene.materials[1]; }
  opengl_shape&    glshape() { return scene.shapes[0]; }
  opengl_shape&    glpoints() { return scene.shapes[1]; }
  opengl_shape&    glvector_field() { return scene.shapes[2]; }
  opengl_shape&    gledges() { return scene.shapes[3]; }
  opengl_shape&    glpolyline() { return scene.shapes[4]; }
};

void init_camera(app_state& app) {
  app.camera = yocto_camera{};
  auto from  = vec3f{0, 0.5, 1.5};
  auto to    = vec3f{0, 0, 0};
  auto up    = vec3f{0, 1, 0};
  // camera.frame     = lookat_frame(from, to, up);
  // camera.yfov      = radians(45);

  app.camera.lens         = 0.02f;
  app.camera.orthographic = false;
  app.camera.aperture     = 0;
  app.camera.frame        = lookat_frame(from, to, up);
  app.camera.film         = {0.036f, 0.015f};
  app.camera.focus        = length(to - from);
  app.camera_focus        = app.camera.focus;
}

void update_glpolyline(app_state& app, const vector<vec3f>& vertices) {
  auto& glshape = app.glpolyline();
  delete_glarraybuffer(glshape.positions);
  delete_glelementbuffer(glshape.lines);
  if (vertices.size()) {
    auto elements = vector<vec2i>(vertices.size() - 1);
    for (int i = 0; i < elements.size(); i++) elements[i] = {i, i + 1};
    init_glarraybuffer(glshape.positions, vertices, false);
    init_glelementbuffer(glshape.lines, elements, false);
  }
}

void update_glpoints(app_state& app, const vector<vec3f>& points) {
  auto& glshape = app.glpoints();
  // delete_glarraybuffer(glshape.positions);
  // delete_glelementbuffer(glshape.points);
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
  delete_glarraybuffer(glshape.positions);
  delete_glelementbuffer(glshape.lines);
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
  delete_glarraybuffer(app.gledges().positions);
  delete_glelementbuffer(app.gledges().lines);

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

vector<vec3f> get_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions) {
  if (path.vertices.empty()) return {};

  auto positions = vector<vec3f>();
  positions.reserve(path.vertices.size() + 1);
  positions.push_back(mesh_positions[path.start]);

  for (int i = 0; i < path.vertices.size() - 1; ++i) {
    auto [edge, face, x] = path.vertices[i];
    auto p0              = mesh_positions[edge.x];
    auto p1              = mesh_positions[edge.y];
    auto position        = (1 - x) * p0 + x * p1;
    positions.push_back(position);
  }

  if (path.end != -1) {
    positions.push_back(mesh_positions[path.end]);
  }
  return positions;
}

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

// This function must be called every time the geometry is updated
void update_glshape(app_state& app) {
  auto& glshape = app.glshape();
  auto& shape   = app.shape;
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

void update_glcamera(opengl_camera& glcamera, const yocto_camera& camera) {
  glcamera.frame  = camera.frame;
  glcamera.yfov   = camera_yfov(camera);
  glcamera.asepct = camera_aspect(camera);
  glcamera.near   = 0.001f;
  glcamera.far    = 10000;
}

// This function must be called every time the geometry is updated
void update_bvh(app_state& app) {
  make_triangles_bvh(app.bvh, app.shape.triangles, app.shape.positions,
      app.shape.radius, false, false);
}

void init_opengl_scene(app_state& app) {
  make_glscene(app.scene);
  update_glcamera(app.scene.cameras.emplace_back(), app.camera);

  // materials[0] is the material of the model
  auto shape_material      = opengl_material{};
  shape_material.diffuse   = {1, 0.2, 0};
  shape_material.roughness = 0.3;
  app.scene.materials.push_back(shape_material);

  // materials[1] is the material for the lines and points
  auto lines_material      = opengl_material{};
  lines_material.emission  = {1, 1, 1};
  lines_material.roughness = 0.0;
  app.scene.materials.push_back(lines_material);

  // shapes[0] is the model
  app.scene.shapes.push_back({});
  update_glshape(app);

  // instances[0] is the model
  app.scene.instances.push_back({identity3x4f, 0, 0, false});

  // shapes[1]/instances[1] are the points
  app.scene.shapes.push_back({});
  app.scene.instances.push_back({identity3x4f, 1, 1, false});

  // shapes[2]/instances[2] is the vector field
  app.scene.shapes.push_back({});
  app.scene.instances.push_back({identity3x4f, 2, 1, false});

  // shapes[3]/instances[3] are mesh edges
  app.scene.shapes.push_back({});
  app.scene.instances.push_back({identity3x4f, -1, 1, false});
  update_gledges(app);

  // shapes[4]/instances[4] is the polyline
  app.scene.shapes.push_back({});
  app.scene.instances.push_back({identity3x4f, 4, 1, false});

  // lights
  app.scene.lights.push_back({{5, 5, 5}, {30, 30, 30}, 0});
  app.scene.lights.push_back({{-5, 5, 5}, {30, 30, 30}, 0});
  app.scene.lights.push_back({{0, 5, -5}, {30, 30, 30}, 0});
}

vec2f get_opengl_mouse_pos_normalized(const opengl_window& win) {
  double mouse_posx, mouse_posy;
  glfwGetCursorPos(win.win, &mouse_posx, &mouse_posy);
  int width, height;
  glfwGetWindowSize(win.win, &width, &height);
  return vec2f{(float)(mouse_posx / width), (float)(mouse_posy / height)};
}

void key_callback(
    GLFWwindow* window, int key, int scancode, int action, int mods) {
  // ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);

  auto  win = (opengl_window*)glfwGetWindowUserPointer(window);
  auto& app = *(app_state*)win->user_ptr;

  bool pressing = action == GLFW_PRESS;

  // ignore release
  if (not pressing) return;

  if (key == GLFW_KEY_ENTER) {
    printf("Enter pressed!\n");
  }

  if (key == GLFW_KEY_ESCAPE) {
    printf("Esc pressed!\n");
    init_camera(app);
    update_glcamera(app.scene.cameras[0], app.camera);
  }

  if (key == GLFW_KEY_RIGHT or key == GLFW_KEY_LEFT) {
    printf("Arrow pressed!\n");
    // arrow keys...
  }

  if (key == GLFW_KEY_A) {
    printf("A pressed!\n");
    // characters...
  }
}

void mouse_button_callback(
    GLFWwindow* window, int button, int action, int mods) {
  auto  win = (opengl_window*)glfwGetWindowUserPointer(window);
  auto& app = *(app_state*)win->user_ptr;

  auto is_key_pressed = [](opengl_window* win, int key) {
    return glfwGetKey(win->win, key) == GLFW_PRESS;
  };

  auto press       = action == GLFW_PRESS;
  auto mouse       = get_opengl_mouse_pos_normalized(*win);
  auto right_click = button == GLFW_MOUSE_BUTTON_RIGHT;
  if (right_click && press) {
    auto  ray = eval_camera(app.camera, mouse, {0.5, 0.5});
    int   face;
    vec2f uv;
    float distance;
    auto  hit = intersect_triangles_bvh(app.bvh, app.shape.triangles,
        app.shape.positions, ray, face, uv, distance);

    if (hit) {
      int vertex = -1;
      if (face >= 0 && face < app.shape.triangles.size()) {
        auto xyz = vec3f{uv.x, uv.y, 1 - uv.x - uv.y};
        int  i   = 0;
        if (xyz.x > xyz.y and xyz.x > xyz.z) i = 1;
        if (xyz.y > xyz.x and xyz.y > xyz.z) i = 2;
        vertex = app.shape.triangles[face][i];
      }
      if (is_key_pressed(win, GLFW_KEY_LEFT_SHIFT)) {
        app.vertex_selection.push_back(vertex);
      } else {
        app.vertex_selection = {vertex};
      }

      printf("clicked vertex: %d\n", vertex);

      auto positions = vector<vec3f>(app.vertex_selection.size());
      for (int i = 0; i < positions.size(); ++i) {
        positions[i] = app.shape.positions[app.vertex_selection[i]];
      }
      update_glpoints(app, positions);
    }
  }
}

void mouse_scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
  auto  win  = (opengl_window*)glfwGetWindowUserPointer(window);
  auto& app  = *(app_state*)win->user_ptr;
  float zoom = yoffset > 0 ? 0.1 : -0.1;
  update_turntable(app.camera.frame, app.camera_focus, zero2f, zoom, zero2f);
  update_glcamera(app.scene.cameras[0], app.camera);
}

void draw_glwidgets(const opengl_window& win, app_state& app) {
  if (draw_glslider(
          win, "vector field scale", app.vector_field_scale, 0, 100)) {
    update_glvector_field(app, app.vector_field, app.vector_field_scale);
  }

  if (draw_glbutton(win, "Geodesic gradient field")) {
    if (app.vertex_selection.size() > 1) {
      app.scalar_field = compute_geodesic_distances(
          app.solver, app.vertex_selection);

      auto compute_gradient = [](const vector<vec3i>&  triangles,
                                  const vector<vec3f>& positions,
                                  const vector<float>& field, int face) {
        auto& t      = triangles[face];
        auto  xy     = positions[t.y] - positions[t.x];
        auto  yz     = positions[t.z] - positions[t.y];
        auto  zx     = positions[t.x] - positions[t.z];
        auto  normal = normalize(cross(zx, xy));
        auto  result = zero3f;
        result += field[t.x] * cross(normal, yz);
        result += field[t.y] * cross(normal, zx);
        result += field[t.z] * cross(normal, xy);
        return result;
      };
      app.vector_field = vector<vec3f>(app.shape.triangles.size());
      for (int i = 0; i < app.shape.triangles.size(); ++i) {
        app.vector_field[i] = compute_gradient(
            app.shape.triangles, app.shape.positions, app.scalar_field, i);
      }
      update_glvector_field(app, app.vector_field, 100);
    }
  }

  if (draw_glbutton(win, "Geodesic distance field")) {
    if (app.vertex_selection.size()) {
      app.scalar_field = compute_geodesic_distances(
          app.solver, app.vertex_selection);
      auto colors = vector<vec4f>(app.scalar_field.size());
      for (int i = 0; i < colors.size(); ++i) {
        colors[i]   = vec4f(app.scalar_field[i]);
        colors[i].w = 1;
      }
      init_glarraybuffer(app.glshape().colors, colors);
    }
  }

  if (draw_glbutton(win, "Compute geodesic path")) {
    if (app.vertex_selection.size() > 1) {
      auto positions = vector<vec3f>();
      auto n         = app.vertex_selection.size();
      for (int i = 0; i < n - (n < 3); ++i) {
        auto from  = app.vertex_selection[i];
        auto to    = app.vertex_selection[(i + 1) % n];
        auto field = compute_geodesic_distances(app.solver, {to});
        for (auto& f : field) f = -f;

        // @Speed: Remove tags from function api to avoid usless allcations.
        auto path = follow_gradient_field(app.shape.triangles,
            app.shape.positions, app.face_adjacency,
            vector<int>(app.shape.triangles.size(), 0), 0, field, from, to);

        positions += get_positions_from_path(path, app.shape.positions);
      }
      update_glpolyline(app, positions);
    }
  }

  if (draw_glcheckbox(win, "Show mesh edges", app.show_edges)) {
    if (app.show_edges)
      app.scene.instances[3].shape = 3;
    else
      app.scene.instances[3].shape = -1;
  }

  if (draw_glbutton(win, "Clear")) {
    for (int i = 1; i < app.scene.shapes.size(); i++) {
      delete_glshape(app.scene.shapes[i]);
    }
    // For points, only delete elements and keep positions.
    app.vertex_selection.clear();
    delete_glarraybuffer(app.glshape().colors);
    init_glarraybuffer(app.glshape().colors,
        vector<vec4f>(app.shape.positions.size(), {1, 1, 1, 1}));
  }
}

// draw with shading
void draw(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  draw_glscene(
      app.scene, get_glframebuffer_viewport(win, false), app.opengl_options);
  begin_glwidgets(win);
  draw_glwidgets(win, app);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void run_app(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimshproc", &app, draw);
  init_opengl_scene(app);

  glfwSetMouseButtonCallback(win.win, mouse_button_callback);
  glfwSetScrollCallback(win.win, mouse_scroll_callback);
  glfwSetKeyCallback(win.win, key_callback);

  // init widget
  init_glwidgets(win);

  // loop
  auto mouse_pos = zero2f, last_pos = zero2f;

  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto alt_down       = get_glalt_key(win);
    auto shift_down     = get_glshift_key(win);
    auto widgets_active = get_glwidgets_active(win);

    // handle mouse and keyboard for navigation
    if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
      // auto& camera = app.scene.cameras.at(app.opengl_options.camera);
      auto& camera = app.camera;
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      // if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
      rotate.y = -rotate.y;
      pan.x    = -pan.x;
      update_turntable(camera.frame, app.camera_focus, rotate, dolly, pan);
      // update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      update_glcamera(app.scene.cameras[0], camera);
    }

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // clear
  delete_glwindow(win);
}

int main(int num_args, const char* args[]) {
  auto app = app_state{};

  // parse command line
  auto cli = make_cli("yimshproc", "interactive viewer for mesh processing");
  add_cli_option(cli, "--resolution,-r", app.opengl_options.resolution,
      "Image resolution.");
  add_cli_option(cli, "scenes", app.input_filename, "Scene filenames", true);
  if (!parse_cli(cli, num_args, args)) exit(1);

  // make default camera
  load_shape(app.input_filename, app.shape.points, app.shape.lines,
      app.shape.triangles, app.shape.quads, app.shape.positions,
      app.shape.normals, app.shape.texcoords, app.shape.colors,
      app.shape.radius);
  app.face_adjacency   = face_adjacencies(app.shape.triangles);
  app.vertex_adjacency = vertex_adjacencies(
      app.shape.triangles, app.face_adjacency);
  app.solver = make_geodesic_solver(
      app.shape.triangles, app.face_adjacency, app.shape.positions);

  update_bvh(app);
  init_camera(app);

  run_app(app);
}
