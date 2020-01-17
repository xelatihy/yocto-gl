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
  function<void(shared_ptr<app_state>, const opengl_window*)>   draw_glwidgets;

  // Geometry data
  sceneio_shape shape;

  // OpenGL data
  unique_ptr<opengl_scene> scene          = {};
  draw_glscene_params      opengl_options = {};

  // Interaction data
  float          time       = 0;
  bool           show_edges = false;
  sceneio_camera camera;
  float          camera_focus;
  bvh_tree       bvh;

  // Internal handles
  int glshape_id, glpoints_id, glvector_field_id, gledges_id, glpolyline_id;
};

void update_glshape(shared_ptr<app_state> app) {
  // @Issue: This app is specialized for a model that is a triangle mesh.
  //    Loading a generic shape is unsafe, maybe we should load only
  //    triangle meshes here...
  auto& shape = app->shape;
  if (!shape.points.empty()) {
    set_shape(app->scene.get(), app->glshape_id, shape.points, shape.positions,
        shape.normals, shape.texcoords, shape.colors);
  } else if (!shape.lines.empty()) {
    set_shape(app->scene.get(), app->glshape_id, shape.lines, shape.positions,
        shape.normals, shape.texcoords, shape.colors);
  } else if (!shape.triangles.empty()) {
    set_shape(app->scene.get(), app->glshape_id, shape.triangles,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.tangents);
  } else if (!shape.quads.empty()) {
    set_shape(app->scene.get(), app->glshape_id, shape.quads, shape.positions,
        shape.normals, shape.texcoords, shape.colors, shape.tangents);
  }
}

void update_glpolyline(
    shared_ptr<app_state> app, const vector<vec3f>& vertices) {
  if (vertices.size()) {
    auto elements = vector<vec2i>(vertices.size() - 1);
    for (int i = 0; i < elements.size(); i++) elements[i] = {i, i + 1};
    set_shape(app->scene.get(), app->glpolyline_id, elements, vertices, {}, {});
  }
}

void update_glpoints(shared_ptr<app_state> app, const vector<vec3f>& points) {
  if (points.size()) {
    auto elements = vector<int>(points.size());
    for (int i = 0; i < elements.size(); i++) elements[i] = i;
    auto normals = vector<vec3f>(points.size(), {0, 0, 1});
    set_shape(
        app->scene.get(), app->glpoints_id, elements, points, normals, {});
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

  set_shape(
      app->scene.get(), app->glvector_field_id, elements, positions, {}, {});
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
  set_shape(app->scene.get(), app->gledges_id, elements, positions, {}, {});
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
  set_instance(app->scene.get(), app->gledges_id, identity3x4f, -1, 1);
}
void show_edges(shared_ptr<app_state> app) {
  app->show_edges = true;
  set_instance(
      app->scene.get(), app->gledges_id, identity3x4f, app->gledges_id, 1);
}

void init_opengl_scene(shared_ptr<app_state> app) {
  app->scene = unique_ptr<opengl_scene>{make_glscene()};
  add_camera(app->scene.get(), app->camera.frame, app->camera.lens,
      app->camera.aspect, app->camera.film, 0.001, 10000);

  auto shape_material = add_material(app->scene.get());
  set_material_diffuse(app->scene.get(), shape_material, {1, 0.2, 0});
  set_material_roughness(app->scene.get(), shape_material, 0.3);

  // @Issue: Right now we're missing APIs to color things easily.
  auto lines_material = add_material(app->scene.get());
  set_material_emission(app->scene.get(), lines_material, {1, 1, 1});
  set_material_roughness(app->scene.get(), lines_material, 0.0);

  // The model.
  app->glshape_id = add_shape(app->scene.get());
  update_glshape(app);

  // The points.
  app->glpoints_id = add_shape(app->scene.get());

  // The vector field.
  app->glvector_field_id = add_shape(app->scene.get());

  // The edges.
  app->gledges_id = add_shape(app->scene.get());
  update_gledges(app);

  // The polyline.
  app->glpolyline_id = add_shape(app->scene.get());

  // Add instances.
  for (int i = 0; i < 5; ++i) {
    add_instance(app->scene.get(), identity3x4f, i, i ? 1 : 0);
  }

  // Hide edges.
  if (!app->show_edges) hide_edges(app);

  // Add lights.
  add_light(app->scene.get(), {5, 5, 5}, {30, 30, 30}, false);
  add_light(app->scene.get(), {-5, 5, 5}, {30, 30, 30}, false);
  add_light(app->scene.get(), {0, 5, -5}, {30, 30, 30}, false);
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
    function<void(shared_ptr<app_state>, const opengl_window* win)>
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
  auto win = make_glwindow({1280 + 320, 720}, "yimshproc", true);
  init_opengl_scene(app);

  // callbacks
  set_draw_glcallback(
      win, [app](const opengl_window* win, const opengl_input& input) {
        draw_glscene(
            app->scene.get(), input.framebuffer_viewport, app->opengl_options);
      });
  set_widgets_glcallback(
      win, [app, draw_glwidgets](const opengl_window* win,
               const opengl_input& input) { draw_glwidgets(app, win); });
  set_click_glcallback(win, [app](const opengl_window* win, bool left,
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
  set_scroll_glcallback(win, [app](const opengl_window* win, float yoffset,
                                 const opengl_input& input) {
    float zoom = yoffset > 0 ? 0.1 : -0.1;
    update_turntable(
        app->camera.frame, app->camera.focus, zero2f, zoom, zero2f);
    set_camera(app->scene.get(), 0, app->camera.frame, app->camera.lens,
        app->camera.aspect, app->camera.film, 0.001, 10000);
  });
  set_key_glcallback(win, [app](const opengl_window* win, int key,
                              bool pressing, const opengl_input& input) {
    app->key_callback(app, key, pressing);
  });
  set_uiupdate_glcallback(
      win, [app](const opengl_window* win, const opengl_input& input) {
        // Handle mouse and keyboard for navigation.
        if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
            !input.widgets_active) {
          auto& camera = app->camera;
          auto  dolly  = 0.0f;
          auto  pan    = zero2f;
          auto  rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) / 100.0f;
          rotate.y = -rotate.y;
          pan.x    = -pan.x;
          update_turntable(camera.frame, app->camera.focus, rotate, dolly, pan);
          set_camera(app->scene.get(), 0, camera.frame, camera.lens,
              camera.aspect, camera.film, 0.001, 10000);
        }
      });

  // cleanup
  delete_glwindow(win);
}
