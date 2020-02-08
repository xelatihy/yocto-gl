#include <memory>

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;
using namespace std;

struct app_state {
  // Callbacks available for user to build its own behaviors
  function<void(shared_ptr<app_state>)>                         init;
  function<void(shared_ptr<app_state>, int, bool)>              key_callback;
  function<void(shared_ptr<app_state>, int, vec2f, int, float)> click_callback;
  function<void(shared_ptr<app_state>, const opengl_window&)>   draw_glwidgets;

  // Geometry data
  sceneio_shape shape;

  // OpenGL data
  opengl_scene        scene          = {};
  draw_glscene_params opengl_options = {};

  // Interaction data
  float          time       = 0;
  bool           show_edges = false;
  sceneio_camera camera;
  float          camera_focus;
  bvh_tree       bvh;

  // Internal handles
  int camera_id = 0;
  int glshape_id, glpoints_id, glvector_field_id, gledges_id, glpolyline_id;
};

void update_glshape(shared_ptr<app_state> app) {
  // @Issue: This app is specialized for a model that is a triangle mesh.
  //    Loading a generic shape is unsafe, maybe we should load only
  //    triangle meshes here...
  auto& shape = app->shape;
  set_shape_points(app->scene, app->glshape_id, shape.points);
  set_shape_lines(app->scene, app->glshape_id, shape.lines);
  set_shape_triangles(app->scene, app->glshape_id, shape.triangles);
  set_shape_quads(app->scene, app->glshape_id, shape.quads);
  set_shape_positions(app->scene, app->glshape_id, shape.positions);
  set_shape_normals(app->scene, app->glshape_id, shape.normals);
  set_shape_texcoords(app->scene, app->glshape_id, shape.texcoords);
  set_shape_colors(app->scene, app->glshape_id, shape.colors);
}

void update_glpolyline(
    shared_ptr<app_state> app, const vector<vec3f>& vertices) {
  if (vertices.size()) {
    auto elements = vector<vec2i>(vertices.size() - 1);
    for (int i = 0; i < elements.size(); i++) elements[i] = {i, i + 1};
    set_shape_positions(app->scene, app->glpolyline_id, vertices);
    set_shape_lines(app->scene, app->glpolyline_id, elements);
  }
}

void update_glpoints(shared_ptr<app_state> app, const vector<vec3f>& points) {
  if (points.size()) {
    auto elements = vector<int>(points.size());
    for (int i = 0; i < elements.size(); i++) elements[i] = i;
    auto normals = vector<vec3f>(points.size(), {0, 0, 1});
    set_shape_positions(app->scene, app->glpolyline_id, points);
    set_shape_points(app->scene, app->glpolyline_id, elements);
  }
}

void update_glvector_field(shared_ptr<app_state> app,
    const vector<vec3f>& vector_field, float scale = 0.01) {
  auto perface   = vector_field.size() == app->shape.triangles.size();
  auto pervertex = vector_field.size() == app->shape.positions.size();

  if (!perface && !pervertex) {
    throw runtime_error("input vector field has wrong size\n");
  }

  auto size = perface ? app->shape.triangles.size()
                      : app->shape.positions.size();
  auto positions = vector<vec3f>(size * 2);

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

  auto elements = vector<vec2i>(size);
  for (int i = 0; i < elements.size(); i++) {
    elements[i] = {2 * i, 2 * i + 1};
  }

  set_shape_positions(app->scene, app->glvector_field_id, positions);
  set_shape_lines(app->scene, app->glvector_field_id, elements);
}

void update_gledges(shared_ptr<app_state> app) {
  auto positions = app->shape.positions;
  for (int i = 0; i < positions.size(); i++) {
    positions[i] += app->shape.normals[i] * 0.0001;
  }

  auto elements = vector<vec2i>();
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
  set_shape_positions(app->scene, app->gledges_id, positions);
  set_shape_lines(app->scene, app->gledges_id, elements);
}

void init_camera(shared_ptr<app_state> app,
    const vec3f& from = vec3f{0, 0.5, 1.5}, const vec3f& to = {0, 0, 0}) {
  app->camera              = sceneio_camera{};
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

void init_bvh(shared_ptr<app_state> app) {
  make_triangles_bvh(
      app->bvh, app->shape.triangles, app->shape.positions, app->shape.radius);
}

void hide_edges(shared_ptr<app_state> app) {
  app->show_edges = false;
  set_shape_hidden(app->scene, app->gledges_id, true);
}
void show_edges(shared_ptr<app_state> app) {
  app->show_edges = true;
  set_shape_hidden(app->scene, app->gledges_id, false);
}

void init_opengl_scene(shared_ptr<app_state> app) {
  init_glscene(app->scene);
  app->camera_id = add_camera(app->scene);
  set_camera_frame(app->scene, app->camera_id, app->camera.frame);
  set_camera_lens(app->scene, app->camera_id, app->camera.lens,
      app->camera.aspect, app->camera.film);
  set_camera_nearfar(app->scene, app->camera_id, 0.001, 10000);

  // The model.
  app->glshape_id = add_shape(app->scene);
  set_shape_color(app->scene, app->glshape_id, {1, 0.2, 0});
  set_shape_specular(app->scene, app->glshape_id, 1);
  set_shape_roughness(app->scene, app->glshape_id, 0.3);
  update_glshape(app);

  // The points.
  app->glpoints_id = add_shape(app->scene);
  set_shape_emission(app->scene, app->glpoints_id, {1, 1, 1});
  set_shape_roughness(app->scene, app->glpoints_id, 0.0);

  // The vector field.
  app->glvector_field_id = add_shape(app->scene);
  set_shape_emission(app->scene, app->glvector_field_id, {1, 1, 1});
  set_shape_roughness(app->scene, app->glvector_field_id, 0.0);

  // The edges.
  app->gledges_id = add_shape(app->scene);
  set_shape_emission(app->scene, app->gledges_id, {1, 1, 1});
  set_shape_roughness(app->scene, app->gledges_id, 0.0);
  update_gledges(app);

  // The polyline.
  app->glpolyline_id = add_shape(app->scene);
  set_shape_emission(app->scene, app->glpolyline_id, {1, 1, 1});
  set_shape_roughness(app->scene, app->glpolyline_id, 0.0);

  // Hide edges.
  if (!app->show_edges) hide_edges(app);

  // Add lights.
  auto l0 = add_light(app->scene);
  set_light(app->scene, l0, {5, 5, 5}, {30, 30, 30}, false);
  auto l1 = add_light(app->scene);
  set_light(app->scene, l1, {-5, 5, 5}, {30, 30, 30}, false);
  auto l2 = add_light(app->scene);
  set_light(app->scene, l2, {0, 5, -5}, {30, 30, 30}, false);
}

void clear(shared_ptr<app_state> app) {
  // TODO: not sure how this works
  // TODO: fix me
  // for (int i = 0; i < app->scene.shapes.size(); i++) {
  //   if (i == app->glshape_id) continue;
  //   if (i == app->gledges_id) continue;
  //   delete_glshape(app->scene.shapes[i]);
  // }
  // delete_glarraybuffer(app->glshape().colors);
  // init_glarraybuffer(app->glshape().colors,
  //     vector<vec4f>(app->shape.positions.size(), {1, 1, 1, 1}));
}

void yimshproc(const string&                         input_filename,
    function<void(shared_ptr<app_state>)>            init,
    function<void(shared_ptr<app_state>, int, bool)> key_callback,
    function<void(shared_ptr<app_state>, int, vec2f, int, float)>
        click_callback,
    function<void(shared_ptr<app_state>, const opengl_window& win)>
        draw_glwidgets) {
  auto app = make_shared<app_state>();

  // init shape
  load_shape(input_filename, app->shape.points, app->shape.lines,
      app->shape.triangles, app->shape.quads, app->shape.positions,
      app->shape.normals, app->shape.texcoords, app->shape.colors,
      app->shape.radius);
  init_bvh(app);
  init_camera(app);

  app->init           = init;
  app->key_callback   = key_callback;
  app->click_callback = click_callback;
  app->draw_glwidgets = draw_glwidgets;  // @Issue: win not needed for widgets

  app->init(app);

  // Init window.
  auto win = opengl_window{};
  init_glwindow(win, {1280 + 320, 720}, "yimshproc", true);
  init_opengl_scene(app);

  // callbacks
  set_draw_glcallback(win, [app](const opengl_window& win,
                               const opengl_input&    input) {
    draw_glscene(app->scene, input.framebuffer_viewport, app->opengl_options);
  });
  set_widgets_glcallback(
      win, [app, draw_glwidgets](const opengl_window& win,
               const opengl_input& input) { draw_glwidgets(app, win); });
  set_click_glcallback(win, [app](const opengl_window& win, bool left,
                                bool press, const opengl_input& input) {
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
  set_scroll_glcallback(win, [app](const opengl_window& win, float yoffset,
                                 const opengl_input& input) {
    float zoom = yoffset > 0 ? 0.1 : -0.1;
    update_turntable(
        app->camera.frame, app->camera.focus, zero2f, zoom, zero2f);
    set_camera_frame(app->scene, app->camera_id, app->camera.frame);
    set_camera_lens(app->scene, app->camera_id, app->camera.lens,
        app->camera.aspect, app->camera.film);
  });
  set_key_glcallback(win, [app](const opengl_window& win, int key,
                              bool pressing, const opengl_input& input) {
    app->key_callback(app, key, pressing);
  });
  set_uiupdate_glcallback(
      win, [app](const opengl_window& win, const opengl_input& input) {
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
          set_camera_frame(app->scene, app->camera_id, app->camera.frame);
          set_camera_lens(app->scene, app->camera_id, app->camera.lens,
              app->camera.aspect, app->camera.film);
        }
      });

  // cleanup
  clear_glwindow(win);
}
