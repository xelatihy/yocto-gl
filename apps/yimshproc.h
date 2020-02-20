#include <memory>

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace ym;

using std::function;
using std::string;
using std::vector;
using namespace std::string_literals;

struct app_state {
  // Callbacks available for user to build its own behaviors
  function<void(app_state*)>                         init;
  function<void(app_state*, int, bool)>              key_callback;
  function<void(app_state*, int, vec2f, int, float)> click_callback;
  function<void(app_state*, yglu::window*)>          draw_widgets;

  // Geometry data
  yscn::shape shape;

  // OpenGL data
  yglu::scene*       glscene        = new yglu::scene{};
  yglu::scene_params opengl_options = {};

  // Interaction data
  float          time       = 0;
  bool           show_edges = false;
  yscn::camera   camera;
  float          camera_focus;
  ybvh::bvh_tree bvh;

  // Internal handles
  yglu::camera*   glcamera    = nullptr;
  yglu::shape*    glshapes    = nullptr;
  yglu::shape*    glpoints    = nullptr;
  yglu::shape*    glvfields   = nullptr;
  yglu::shape*    gledges     = nullptr;
  yglu::shape*    glpolylines = nullptr;
  yglu::material* glshapem    = nullptr;
  yglu::material* glpointm    = nullptr;
  yglu::material* glvfieldm   = nullptr;
  yglu::material* gledgem     = nullptr;
  yglu::material* glpolylinem = nullptr;
  yglu::object*   glshapeo    = nullptr;
  yglu::object*   glpointo    = nullptr;
  yglu::object*   glvfieldo   = nullptr;
  yglu::object*   gledgeo     = nullptr;
  yglu::object*   glpolylineo = nullptr;

  // cleanup
  ~app_state() {
    if (glscene) delete glscene;
  }
};

void update_glshape(app_state* app) {
  // @Issue: This app is specialized for a model that is a triangle mesh.
  //    Loading a generic shape is unsafe, maybe we should load only
  //    triangle meshes here...
  auto& shape = app->shape;
  set_points(app->glshapes, shape.points);
  set_lines(app->glshapes, shape.lines);
  set_triangles(app->glshapes, shape.triangles);
  set_quads(app->glshapes, shape.quads);
  set_positions(app->glshapes, shape.positions);
  set_normals(app->glshapes, shape.normals);
  set_texcoords(app->glshapes, shape.texcoords);
  set_colors(app->glshapes, shape.colors);
}

void update_glpolyline(app_state* app, const std::vector<vec3f>& vertices) {
  if (vertices.size()) {
    auto elements = std::vector<vec2i>(vertices.size() - 1);
    for (int i = 0; i < elements.size(); i++) elements[i] = {i, i + 1};
    set_positions(app->glpolylines, vertices);
    set_lines(app->glpolylines, elements);
  }
}

void update_glpoints(app_state* app, const std::vector<vec3f>& points) {
  if (points.size()) {
    auto elements = std::vector<int>(points.size());
    for (int i = 0; i < elements.size(); i++) elements[i] = i;
    auto normals = std::vector<vec3f>(points.size(), {0, 0, 1});
    set_positions(app->glpolylines, points);
    set_points(app->glpolylines, elements);
  }
}

void update_glvector_field(app_state* app,
    const std::vector<vec3f>& vector_field, float scale = 0.01) {
  auto perface   = vector_field.size() == app->shape.triangles.size();
  auto pervertex = vector_field.size() == app->shape.positions.size();

  if (!perface && !pervertex) {
    throw std::runtime_error("input vector field has wrong size\n");
  }

  auto size = perface ? app->shape.triangles.size()
                      : app->shape.positions.size();
  auto positions = std::vector<vec3f>(size * 2);

  // Per-face vector field
  if (perface) {
    for (int i = 0; i < app->shape.triangles.size(); i++) {
      auto x      = app->shape.positions[app->shape.triangles[i].x];
      auto y      = app->shape.positions[app->shape.triangles[i].y];
      auto z      = app->shape.positions[app->shape.triangles[i].z];
      auto normal = triangle_normal(x, y, z);
      normal *= 0.0001;
      auto center          = (x + y + z) / 3;
      auto from            = center + normal;
      auto to              = from + (scale * vector_field[i]) + normal;
      positions[i * 2]     = from;
      positions[i * 2 + 1] = to;
    }
  } else {
    for (int i = 0; i < app->shape.positions.size(); i++) {
      auto from            = app->shape.positions[i];
      auto to              = from + scale * vector_field[i];
      positions[i * 2]     = from;
      positions[i * 2 + 1] = to;
    }
  }

  auto elements = std::vector<vec2i>(size);
  for (int i = 0; i < elements.size(); i++) {
    elements[i] = {2 * i, 2 * i + 1};
  }

  set_positions(app->glvfields, positions);
  set_lines(app->glvfields, elements);
}

void update_gledges(app_state* app) {
  auto positions = app->shape.positions;
  for (int i = 0; i < positions.size(); i++) {
    positions[i] += app->shape.normals[i] * 0.0001;
  }

  auto elements = std::vector<vec2i>();
  elements.reserve(app->shape.triangles.size() * 3);
  for (int i = 0; i < app->shape.triangles.size(); i++) {
    for (int k = 0; k < 3; k++) {
      auto a = app->shape.triangles[i][k];
      auto b = app->shape.triangles[i][(k + 1) % 3];
      if (a < b) {
        elements.push_back({a, b});
      }
    }
  }
  set_positions(app->gledges, positions);
  set_lines(app->gledges, elements);
}

void init_camera(app_state* app, const vec3f& from = vec3f{0, 0.5, 1.5},
    const vec3f& to = {0, 0, 0}) {
  app->camera              = yscn::camera{};
  auto up                  = vec3f{0, 1, 0};
  app->camera.lens         = 0.02f;
  app->camera.orthographic = false;
  app->camera.aperture     = 0;
  app->camera.frame        = lookat_frame(from, to, up);
  app->camera.film         = 0.036f;
  app->camera.aspect       = 0.036f / 0.015f;
  app->camera.focus        = length(to - from);
  app->camera_focus        = app->camera.focus;
}

void init_bvh(app_state* app) {
  make_triangles_bvh(
      app->bvh, app->shape.triangles, app->shape.positions, app->shape.radius);
}

void hide_edges(app_state* app) {
  app->show_edges = false;
  set_hidden(app->gledgeo, true);
}
void show_edges(app_state* app) {
  app->show_edges = true;
  set_hidden(app->gledgeo, false);
}

void init_opengl_scene(app_state* app) {
  init_glscene(app->glscene);
  app->glcamera = add_camera(app->glscene);
  set_frame(app->glcamera, app->camera.frame);
  set_lens(
      app->glcamera, app->camera.lens, app->camera.aspect, app->camera.film);
  set_nearfar(app->glcamera, 0.001, 10000);

  // The model.
  app->glshapes = add_shape(app->glscene);
  app->glshapem = add_material(app->glscene);
  set_color(app->glshapem, {1, 0.2, 0});
  set_specular(app->glshapem, 1);
  set_roughness(app->glshapem, 0.3);
  app->glshapeo = add_object(app->glscene);
  set_shape(app->glshapeo, app->glshapes);
  set_material(app->glshapeo, app->glshapem);
  update_glshape(app);

  // The points.
  app->glpoints = add_shape(app->glscene);
  app->glpointm = add_material(app->glscene);
  set_emission(app->glpointm, {1, 1, 1});
  set_roughness(app->glpointm, 0.0);
  app->glpointo = add_object(app->glscene);
  set_shape(app->glpointo, app->glpoints);
  set_material(app->glpointo, app->glpointm);

  // The vector field.
  app->glvfields = add_shape(app->glscene);
  app->glvfieldm = add_material(app->glscene);
  set_emission(app->glvfieldm, {1, 1, 1});
  set_roughness(app->glvfieldm, 0.0);
  app->glvfieldo = add_object(app->glscene);
  set_shape(app->glvfieldo, app->glvfields);
  set_material(app->glvfieldo, app->glvfieldm);

  // The edges.
  app->gledges = add_shape(app->glscene);
  app->gledgem = add_material(app->glscene);
  set_emission(app->gledgem, {1, 1, 1});
  set_roughness(app->gledgem, 0.0);
  app->gledgeo = add_object(app->glscene);
  set_shape(app->gledgeo, app->gledges);
  set_material(app->gledgeo, app->glvfieldm);
  update_gledges(app);

  // The polyline.
  app->glpolylines = add_shape(app->glscene);
  app->glpolylinem = add_material(app->glscene);
  set_emission(app->glpolylinem, {1, 1, 1});
  set_roughness(app->glpolylinem, 0.0);
  app->glpolylineo = add_object(app->glscene);
  set_shape(app->glpolylineo, app->glpolylines);
  set_material(app->glpolylineo, app->glpolylinem);

  // Hide edges.
  if (!app->show_edges) hide_edges(app);

  // Add lights.
  set_light(add_light(app->glscene), {5, 5, 5}, {30, 30, 30}, false);
  set_light(add_light(app->glscene), {-5, 5, 5}, {30, 30, 30}, false);
  set_light(add_light(app->glscene), {0, 5, -5}, {30, 30, 30}, false);
}

void clear(app_state* app) {
  // TODO: not sure how this works
  // TODO: fix me
  // for (int i = 0; i < app->scene.shapes.size(); i++) {
  //   if (i == app->glshape_id) continue;
  //   if (i == app->gledges_id) continue;
  //   delete_glshape(app->scene.shapes[i]);
  // }
  // delete_glarraybuffer(app->glshapes().colors);
  // init_glarraybuffer(app->glshapes().colors,
  //     std::vector<vec4f>(app->shape.positions.size(), {1, 1, 1, 1}));
}

void yimshproc(const std::string&                      input_filename,
    function<void(app_state*)>                         init,
    function<void(app_state*, int, bool)>              key_callback,
    function<void(app_state*, int, vec2f, int, float)> click_callback,
    function<void(app_state*, yglu::window* win)>      draw_widgets) {
  auto app_guard = std::make_unique<app_state>();
  auto app       = app_guard.get();

  // init shape
  auto ioerror = ""s;
  if (!yshp::load_shape(input_filename, app->shape.points, app->shape.lines,
          app->shape.triangles, app->shape.quads, app->shape.positions,
          app->shape.normals, app->shape.texcoords, app->shape.colors,
          app->shape.radius, ioerror))
    ycli::print_fatal(ioerror);
  init_bvh(app);
  init_camera(app);

  app->init           = init;
  app->key_callback   = key_callback;
  app->click_callback = click_callback;
  app->draw_widgets   = draw_widgets;  // @Issue: win not needed for widgets

  app->init(app);

  // Init window.
  auto win_guard = std::make_unique<yglu::window>();
  auto win       = win_guard.get();
  init_glwindow(win, {1280 + 320, 720}, "yimshproc", true);
  init_opengl_scene(app);

  // callbacks
  set_draw_callback(win, [app](yglu::window* win, const yglu::input& input) {
    draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
        app->opengl_options);
  });
  set_widgets_callback(
      win, [app, draw_widgets](yglu::window* win, const yglu::input& input) {
        draw_widgets(app, win);
      });
  set_click_callback(win, [app](yglu::window* win, bool left, bool press,
                              const yglu::input& input) {
    auto mouse = input.mouse_pos /
                 vec2f{(float)input.window_size.x, (float)input.window_size.y};

    // Ray trace camera ray.
    if (!left && press) {
      auto ray          = camera_ray(app->camera.frame, app->camera.lens,
          app->camera.aspect >= 1
              ? vec2f{app->camera.film, app->camera.film / app->camera.aspect}
              : vec2f{app->camera.film * app->camera.aspect, app->camera.film},
          mouse + 0.5f);
      auto intersection = intersect_triangles_bvh(
          app->bvh, app->shape.triangles, app->shape.positions, ray);

      if (intersection.hit) {
        auto uvw = vec3f{intersection.uv.x, intersection.uv.y,
            1 - intersection.uv.x - intersection.uv.y};
        int  k   = 0;
        if (uvw.x > uvw.y && uvw.x > uvw.z) k = 1;
        if (uvw.y > uvw.x && uvw.y > uvw.z) k = 2;
        auto vertex = app->shape.triangles[intersection.element][k];
        app->click_callback(app, intersection.element, intersection.uv, vertex,
            intersection.distance);
      }
    }
  });
  set_scroll_callback(
      win, [app](yglu::window* win, float yoffset, const yglu::input& input) {
        float zoom = yoffset > 0 ? 0.1 : -0.1;
        update_turntable(
            app->camera.frame, app->camera.focus, zero2f, zoom, zero2f);
        set_frame(app->glcamera, app->camera.frame);
        set_lens(app->glcamera, app->camera.lens, app->camera.aspect,
            app->camera.film);
      });
  set_key_callback(win,
      [app](yglu::window* win, int key, bool pressing,
          const yglu::input& input) { app->key_callback(app, key, pressing); });
  set_uiupdate_callback(
      win, [app](yglu::window* win, const yglu::input& input) {
        // Handle mouse and keyboard for navigation.
        if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
            !input.widgets_active) {
          auto dolly  = 0.0f;
          auto pan    = zero2f;
          auto rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) / 100.0f;
          rotate.y = -rotate.y;
          pan.x    = -pan.x;
          update_turntable(
              app->camera.frame, app->camera.focus, rotate, dolly, pan);
          set_frame(app->glcamera, app->camera.frame);
          set_lens(app->glcamera, app->camera.lens, app->camera.aspect,
              app->camera.film);
        }
      });

  // cleanup
  clear_glwindow(win);
}
