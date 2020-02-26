#include "yimshproc.h"

struct my_data {
  std::vector<vec3i>       face_adjacency;
  std::vector<vector<int>> vertex_adjacency;
  yshp::geodesic_solver    solver;
  std::vector<float>       scalar_field;
  std::vector<vec3f>       vector_field;

  int              num_samples      = 500;
  std::vector<int> vertex_selection = {};
};

void my_init(my_data& data, app_state* app) {
  data.face_adjacency   = yshp::face_adjacencies(app->shape.triangles);
  data.vertex_adjacency = yshp::vertex_adjacencies(
      app->shape.triangles, data.face_adjacency);
  data.solver = yshp::make_geodesic_solver(
      app->shape.triangles, data.face_adjacency, app->shape.positions);
}

void my_keycallback(my_data& data, app_state* app, int key, bool pressing) {
  // Ignore release.
  if (!pressing) return;

  ycli::print_info(
      "press: " + std::to_string((char)key) + " [" + std::to_string(key) + "]");
  auto enter = 257;
  auto esc   = 256;

  if (key == enter) {
    ycli::print_info("Enter pressed!");
  }

  if (key == esc) {
    ycli::print_info("Esc pressed!");
    init_camera(app);
    set_frame(app->glcamera, app->camera.frame);
    set_lens(
        app->glcamera, app->camera.lens, app->camera.aspect, app->camera.film);
  }

  if (key == 'z') {
    ycli::print_info("Z pressed!");
  }
}

void my_click_callback(my_data& data, app_state* app, int face, const vec2f& uv,
    int vertex, float distance) {
  ycli::print_info("clicked vertex: " + std::to_string(vertex));
  data.vertex_selection.push_back(vertex);

  auto positions = std::vector<vec3f>(data.vertex_selection.size());
  for (int i = 0; i < positions.size(); ++i) {
    positions[i] = app->shape.positions[data.vertex_selection[i]];
  }
  update_glpoints(app, positions);
}

void my_draw_glwidgets(my_data& data, app_state* app, ygui::window* win) {
  if (draw_button(win, "Geodesic gradient field")) {
    if (data.vertex_selection.size() > 1) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, data.vertex_selection);

      data.vector_field = std::vector<vec3f>(app->shape.triangles.size());
      for (int i = 0; i < app->shape.triangles.size(); ++i) {
        data.vector_field[i] = yshp::compute_gradient(
            app->shape.triangles[i], app->shape.positions, data.scalar_field);
      }
      update_glvector_field(app, data.vector_field, 100);
    }
  }

  if (draw_button(win, "Geodesic distance field")) {
    if (data.vertex_selection.size()) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, data.vertex_selection);
      auto colors = std::vector<vec3f>(data.scalar_field.size());
      for (int i = 0; i < colors.size(); ++i) {
        colors[i] = vec3f(data.scalar_field[i]);
      }
      set_colors(app->glshapes, colors);
    }
  }

  if (draw_button(win, "Compute geodesic path")) {
    if (data.vertex_selection.size() > 1) {
      auto positions = std::vector<vec3f>();
      auto n         = data.vertex_selection.size();
      for (int i = 0; i < n - (n < 3); ++i) {
        auto from  = data.vertex_selection[i];
        auto to    = data.vertex_selection[(i + 1) % n];
        auto field = compute_geodesic_distances(data.solver, {to});
        for (auto& f : field) f = -f;

        // @Speed: Remove tags from function api to avoid this.
        auto dummy_tags = std::vector<int>(app->shape.triangles.size(), 0);
        auto path       = yshp::integrate_field(app->shape.triangles,
            app->shape.positions, data.face_adjacency, dummy_tags, 0, field,
            from, to);

        auto ppositions = yshp::make_positions_from_path(
            path, app->shape.positions);
        positions.insert(positions.end(), ppositions.begin(), ppositions.end());
      }
      update_glpolyline(app, positions);
    }
  }

  if (draw_slider(win, "Sample points", data.num_samples, 0, 10000)) {
    data.vertex_selection = sample_vertices_poisson(
        data.solver, data.num_samples);
    auto positions = std::vector<vec3f>(data.vertex_selection.size());
    for (int i = 0; i < positions.size(); ++i)
      positions[i] = app->shape.positions[data.vertex_selection[i]];
    update_glpoints(app, positions);
  }

  if (draw_checkbox(win, "Show mesh edges", app->show_edges)) {
    if (app->show_edges)
      show_edges(app);
    else
      hide_edges(app);
  }

  if (draw_button(win, "Clear")) {
    data.vertex_selection.clear();
    clear(app);
  }
}

int main(int argc, const char* argv[]) {
  std::string input_filename = "model.obj";

  // Parse command line.
  auto cli = ycli::make_cli(
      "yimshproc", "interactive viewer for mesh processing");
  add_option(cli, "model", input_filename, "model filenames", true);
  parse_cli(cli, argc, argv);

  auto data = my_data{};

  // Create callbacks that interface with yimshproc.
  auto init = [&data](app_state* app) {
    ycli::print_info("init my data");
    my_init(data, app);
  };
  auto key_callback = [&data](app_state* app, int key, bool pressing) {
    my_keycallback(data, app, key, pressing);
  };
  auto click_callback = [&data](app_state* a, int f, vec2f uv, int v, float d) {
    my_click_callback(data, a, f, uv, v, d);
  };
  auto draw_widgets = [&data](app_state* app, ygui::window* win) {
    my_draw_glwidgets(data, app, win);
  };

  yimshproc(input_filename, init, key_callback, click_callback, draw_widgets);

  // done
  return 0;
}
