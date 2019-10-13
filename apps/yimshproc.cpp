#include "yimshproc.h"
using namespace yocto;

struct my_data {
  vector<vec3i>       face_adjacency;
  vector<vector<int>> vertex_adjacency;
  geodesic_solver     solver;
  vector<float>       scalar_field;
  vector<vec3f>       vector_field;
  int                 num_samples = 500;
};

void my_init(my_data& data, app_state& app) {
  data.face_adjacency   = face_adjacencies(app.shape.triangles);
  data.vertex_adjacency = vertex_adjacencies(
      app.shape.triangles, data.face_adjacency);
  data.solver = make_geodesic_solver(
      app.shape.triangles, data.face_adjacency, app.shape.positions);
}

void my_keycallback(my_data& data, app_state& app, int key, int scancode,
    int action, int mods) {
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

  if (key == GLFW_KEY_A) {
    printf("A pressed!\n");
    // characters...
  }
}

void my_draw_glwidgets(
    my_data& data, app_state& app, const opengl_window& win) {
  if (draw_glbutton(win, "Geodesic gradient field")) {
    if (app.vertex_selection.size() > 1) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, app.vertex_selection);

      // We should make this function public from yocto_shape
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
      data.vector_field = vector<vec3f>(app.shape.triangles.size());
      for (int i = 0; i < app.shape.triangles.size(); ++i) {
        data.vector_field[i] = compute_gradient(
            app.shape.triangles, app.shape.positions, data.scalar_field, i);
      }
      update_glvector_field(app, data.vector_field, 100);
    }
  }

  if (draw_glbutton(win, "Geodesic distance field")) {
    if (app.vertex_selection.size()) {
      data.scalar_field = compute_geodesic_distances(
          data.solver, app.vertex_selection);
      auto colors = vector<vec4f>(data.scalar_field.size());
      for (int i = 0; i < colors.size(); ++i) {
        colors[i]   = vec4f(data.scalar_field[i]);
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
        auto field = compute_geodesic_distances(data.solver, {to});
        for (auto& f : field) f = -f;

        // @Speed: Remove tags from function api to avoid useless allcations.
        auto path = follow_gradient_field(app.shape.triangles,
            app.shape.positions, data.face_adjacency,
            vector<int>(app.shape.triangles.size(), 0), 0, field, from, to);

        positions += make_positions_from_path(path, app.shape.positions);
      }
      update_glpolyline(app, positions);
    }
  }

  if (draw_glslider(win, "Sample points", data.num_samples, 0, 10000)) {
    app.vertex_selection = sample_vertices_poisson(
        data.solver, data.num_samples);
    auto positions = vector<vec3f>(app.vertex_selection.size());
    for (int i = 0; i < positions.size(); ++i)
      positions[i] = app.shape.positions[app.vertex_selection[i]];
    update_glpoints(app, positions);
  }

  if (draw_glcheckbox(win, "Show mesh edges", app.show_edges)) {
    if (app.show_edges)
      show_edges(app);
    else
      hide_edges(app);
  }

  if (draw_glbutton(win, "Clear")) {
    clear(app);
  }
}

int main(int num_args, const char* args[]) {
  string input_filename = "model.obj";

  // parse command line
  auto cli = make_cli("yimshproc", "interactive viewer for mesh processing");
  add_cli_option(cli, "Model", input_filename, "Model filenames", true);
  if (!parse_cli(cli, num_args, args)) exit(1);

  auto data = my_data{};

  // Create callbacks that interface with yimshproc
  auto init         = [&data](app_state& app) { my_init(data, app); };
  auto key_callback = [&data](app_state& app, int key, int s, int a, int m) {
    my_keycallback(data, app, key, s, a, m);
  };
  auto draw_glwidgets = [&data](app_state& app, const opengl_window& win) {
    my_draw_glwidgets(data, app, win);
  };

  yimshproc(input_filename, init, key_callback, draw_glwidgets);
}
