#include "yimshproc.h"
using namespace yocto;

void my_keycallback(
    app_state& app, int key, int scancode, int action, int mods) {
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

void my_draw_glwidgets(app_state& app, const opengl_window& win) {
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

        // @Speed: Remove tags from function api to avoid useless allcations.
        auto path = follow_gradient_field(app.shape.triangles,
            app.shape.positions, app.face_adjacency,
            vector<int>(app.shape.triangles.size(), 0), 0, field, from, to);

        positions += get_positions_from_path(path, app.shape.positions);
      }
      update_glpolyline(app, positions);
    }
  }

  if (draw_glslider(win, "Sample points", app.num_samples, 0, 10000)) {
    app.vertex_selection = sample_vertices_poisson(app.solver, app.num_samples);
    auto positions       = vector<vec3f>(app.vertex_selection.size());
    for (int i = 0; i < positions.size(); ++i)
      positions[i] = app.shape.positions[app.vertex_selection[i]];
    update_glpoints(app, positions);
  }

  if (draw_glcheckbox(win, "Show mesh edges", app.show_edges)) {
    if (app.show_edges)
      app.scene.instances[app.gledges_id].shape = app.gledges_id;
    else
      app.scene.instances[app.gledges_id].shape = -1;
  }

  if (draw_glbutton(win, "Clear")) {
    for (int i = 0; i < app.scene.shapes.size(); i++) {
      if (i == app.glshape_id) continue;
      if (i == app.gledges_id) continue;
      delete_glshape(app.scene.shapes[i]);
    }

    app.vertex_selection.clear();
    delete_glarraybuffer(app.glshape().colors);
    init_glarraybuffer(app.glshape().colors,
        vector<vec4f>(app.shape.positions.size(), {1, 1, 1, 1}));
  }
}

void my_init(app_state& app) {}

int main(int num_args, const char* args[]) {
  string input_filename = "model.obj";

  // parse command line
  auto cli = make_cli("yimshproc", "interactive viewer for mesh processing");
  add_cli_option(cli, "Model", input_filename, "Model filenames", true);
  if (!parse_cli(cli, num_args, args)) exit(1);

  yimshproc(input_filename, my_init, my_keycallback, my_draw_glwidgets);
}
