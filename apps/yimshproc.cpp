#include "yimshproc.h"
using namespace yocto;

struct my_data {
  vector<vec3i>       face_adjacency;
  vector<vector<int>> vertex_adjacency;
  geodesic_solver     solver;
  vector<float>       scalar_field;
  vector<vec3f>       vector_field;

  int         num_samples      = 500;
  vector<int> vertex_selection = {};
};

void my_init(my_data& data, shared_ptr<app_state> app) {
  data.face_adjacency   = face_adjacencies(app->shape.triangles);
  data.vertex_adjacency = vertex_adjacencies(
      app->shape.triangles, data.face_adjacency);
  data.solver = make_geodesic_solver(
      app->shape.triangles, data.face_adjacency, app->shape.positions);
}

void my_keycallback(
    my_data& data, shared_ptr<app_state> app, int key, bool pressing) {
  // Ignore release.
  if (!pressing) return;

  printf("press: %c [%d]\n", (char)key, key);
  auto enter = 257;
  auto esc   = 256;

  if (key == enter) {
    printf("Enter pressed!\n");
  }

  if (key == esc) {
    printf("Esc pressed!\n");
    init_camera(app);
    set_camera_frame(app->scene, app->camera_id, app->camera.frame);
    set_camera_lens(app->scene, app->camera_id, app->camera.lens);
    set_camera_aspect(app->scene, app->camera_id, app->camera.aspect);
    set_camera_film(app->scene, app->camera_id, app->camera.film);
  }

  if (key == 'z') {
    printf("Z pressed!\n");
  }
}

void my_click_callback(my_data& data, shared_ptr<app_state> app, int face,
    const vec2f& uv, int vertex, float distance) {
  printf("clicked vertex: %d\n", vertex);
  data.vertex_selection.push_back(vertex);

  auto positions = vector<vec3f>(data.vertex_selection.size());
  for (int i = 0; i < positions.size(); ++i) {
    positions[i] = app->shape.positions[data.vertex_selection[i]];
  }
  update_glpoints(app, positions);
}

void my_draw_glwidgets(
    my_data& data, shared_ptr<app_state> app, const opengl_window& win) {
  if (draw_glbutton(win, "Geodesic gradient field")) {
    if (data.vertex_selection.size() > 1) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, data.vertex_selection);

      data.vector_field = vector<vec3f>(app->shape.triangles.size());
      for (int i = 0; i < app->shape.triangles.size(); ++i) {
        data.vector_field[i] = compute_gradient(
            app->shape.triangles[i], app->shape.positions, data.scalar_field);
      }
      update_glvector_field(app, data.vector_field, 100);
    }
  }

  if (draw_glbutton(win, "Geodesic distance field")) {
    if (data.vertex_selection.size()) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, data.vertex_selection);
      auto colors = vector<vec4f>(data.scalar_field.size());
      for (int i = 0; i < colors.size(); ++i) {
        colors[i]   = vec4f(data.scalar_field[i]);
        colors[i].w = 1;
      }
      set_shape_colors(app->scene, app->glshape_id, colors);
    }
  }

  if (draw_glbutton(win, "Compute geodesic path")) {
    if (data.vertex_selection.size() > 1) {
      auto positions = vector<vec3f>();
      auto n         = data.vertex_selection.size();
      for (int i = 0; i < n - (n < 3); ++i) {
        auto from  = data.vertex_selection[i];
        auto to    = data.vertex_selection[(i + 1) % n];
        auto field = compute_geodesic_distances(data.solver, {to});
        for (auto& f : field) f = -f;

        // @Speed: Remove tags from function api to avoid this.
        auto dummy_tags = vector<int>(app->shape.triangles.size(), 0);
        auto path = integrate_field(app->shape.triangles, app->shape.positions,
            data.face_adjacency, dummy_tags, 0, field, from, to);

        auto ppositions = make_positions_from_path(path, app->shape.positions);
        positions.insert(positions.end(), ppositions.begin(), ppositions.end());
      }
      update_glpolyline(app, positions);
    }
  }

  if (draw_glslider(win, "Sample points", data.num_samples, 0, 10000)) {
    data.vertex_selection = sample_vertices_poisson(
        data.solver, data.num_samples);
    auto positions = vector<vec3f>(data.vertex_selection.size());
    for (int i = 0; i < positions.size(); ++i)
      positions[i] = app->shape.positions[data.vertex_selection[i]];
    update_glpoints(app, positions);
  }

  if (draw_glcheckbox(win, "Show mesh edges", app->show_edges)) {
    if (app->show_edges)
      show_edges(app);
    else
      hide_edges(app);
  }

  if (draw_glbutton(win, "Clear")) {
    data.vertex_selection.clear();
    clear(app);
  }
}

void run_app(int argc, const char* argv[]) {
  string input_filename = "model.obj";

  // Parse command line.
  auto cli = make_cli("yimshproc", "interactive viewer for mesh processing");
  add_cli_option(cli, "model", input_filename, "model filenames", true);
  parse_cli(cli, argc, argv);

  auto data = my_data{};

  // Create callbacks that interface with yimshproc.
  auto init = [&data](shared_ptr<app_state> app) {
    auto timer = print_timed("init my data");
    my_init(data, app);
    print_elapsed(timer);
  };
  auto key_callback = [&data](
                          shared_ptr<app_state> app, int key, bool pressing) {
    my_keycallback(data, app, key, pressing);
  };
  auto click_callback = [&data](shared_ptr<app_state> a, int f, vec2f uv, int v,
                            float d) {
    my_click_callback(data, a, f, uv, v, d);
  };
  auto draw_glwidgets = [&data](shared_ptr<app_state> app,
                            const opengl_window&      win) {
    my_draw_glwidgets(data, app, win);
  };

  yimshproc(input_filename, init, key_callback, click_callback, draw_glwidgets);
}

int main(int argc, const char* argv[]) {
  try {
    run_app(argc, argv);
    return 0;
  } catch (std::exception& e) {
    print_fatal(e.what());
    return 1;
  }
}
