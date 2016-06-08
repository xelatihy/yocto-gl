//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

// clang-format off
#ifndef __APPLE__
#include "ext/glew/glew.h"
#else
#include "OpenGL/gl.h"
#endif
#include "ext/glfw/glfw3.h"
// clang-format on

#define YGL_DECLARATION

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

// view data
struct view_params {
    // filenames
    ym_string filename;
    ym_string imfilename;

    // window data
    int w = 0, h = 0;
    float exposure = 0, gamma = 2.2f;
    int background = 0;
    ym_vec4f backgrounds[4] = {{0, 0, 0, 0},
                               {0.18f, 0.18f, 0.18f, 0},
                               {0.5f, 0.5f, 0.5f, 0},
                               {1, 1, 1, 0}};

    // drawing
    bool wireframe = false;
    bool edges = false;

    // mouse state
    int mouse_button = 0;
    ym_vec2f mouse_pos;

    // scene
    yo_scene* scene = nullptr;
    ym_vec3f amb = ym_zero3f;
    int cur_camera = 0;

    // shading
    int shade_prog = 0;
    ym_vector<int> shade_txt;

    // lighting
    bool camera_lights = false;
};

view_params* init_view_params(const ym_string& filename,
                              const ym_string& imfilename, yo_scene* scene,
                              int w, int h, bool camera_lights, float amb) {
    view_params* view = new view_params();

    view->filename = filename;
    view->imfilename = imfilename;

    view->w = w;
    view->h = h;

    view->scene = scene;
    view->amb = {amb, amb, amb};

    view->camera_lights = camera_lights;

    return view;
}

void shade(yo_scene* scene, int cur_camera, int prog, int* txt, float exposure,
           float gamma_, bool wireframe, bool edges, bool camera_lights) {
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    yo_camera* cam = &scene->cameras[cur_camera];
    ym_frame3f camera_frame = ym_lookat_xform3(cam->from, cam->to, cam->up);
    ym_mat4f camera_xform = ym_mat4f(camera_frame);
    ym_mat4f camera_view = ym_mat4f(ym_inverse(camera_frame));
    ym_mat4f camera_proj = ym_perspective_mat4(
        2 * atanf(cam->height / 2), cam->width / cam->height, 0.1f, 10000.0f);

    yg_stdshader_begin_frame(prog, camera_lights, exposure, gamma_,
                             camera_xform.data(), camera_view.data(),
                             camera_proj.data());

    if (!camera_lights) {
        int nlights = 0;
        ym_vec3f light_pos[16], light_ke[16];
        int light_type[16];
        for (int i = 0; i < scene->shapes.size() && nlights < 16; i++) {
            yo_shape* shape = &scene->shapes[i];
            if (shape->etype != yo_etype_point) continue;
            if (shape->matid < 0) continue;
            yo_material* mat = &scene->materials[shape->matid];
            if (!mat->ke[0] && !mat->ke[1] && !mat->ke[2]) continue;
            for (int j = 0; j < shape->nverts && nlights < 16; j++) {
                light_pos[nlights] = shape->pos[j];
                if (shape->xformed) {
                    light_pos[nlights] =
                        ym_transform_point(shape->xform, light_pos[nlights]);
                }
                light_ke[nlights] = mat->ke;
                light_type[nlights] = 0;
                nlights++;
            }
        }
        ym_vec3f zero = {0, 0, 0};
        yg_stdshader_set_lights(prog, &zero.x, nlights, (float*)light_pos,
                                (float*)light_ke, light_type);
    }

    for (int i = 0; i < scene->shapes.size(); i++) {
        yo_shape* shape = &scene->shapes[i];

        yg_stdshader_begin_shape(prog, ym_mat4f(shape->xform).data());

        if (shape->matid >= 0) {
            yo_material* mat = &scene->materials[shape->matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yg_stdshader_set_material(
                prog, shape->etype, &mat->ke.x, &mat->kd.x, &mat->ks.x,
                yg_specular_exponent_to_roughness(mat->ns),
                __txt(mat->ke_txtid), __txt(mat->kd_txtid),
                __txt(mat->ks_txtid), __txt(mat->ns_txtid), false);
        } else {
            ym_vec3f kd = {0.8f, 0.8f, 0.8f};
            ym_vec3f zero = ym_zero3f;
            yg_stdshader_set_material(prog, shape->etype, &zero.x, &kd.x,
                                      &zero.x, 0, 0, 0, 0, 0, false);
        }

        yg_stdshader_set_vert(prog, (float*)shape->pos.data(),
                              (float*)shape->norm.data(),
                              (float*)shape->texcoord.data(), 0);

        yg_stdshader_draw_elem(prog, shape->nelems, shape->elem.data(),
                               shape->etype);

        if (edges && !wireframe) {
            float zero[3] = {0, 0, 0};
            yg_stdshader_set_material(prog, shape->etype, zero, zero, zero, 0,
                                      0, 0, 0, 0, false);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yg_stdshader_draw_elem(prog, shape->nelems, shape->elem.data(),
                                   shape->etype);
            glDepthRange(0, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glLineWidth(1);
        }

        yg_stdshader_end_shape();
    }

    yg_stdshader_end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void save_screenshot(GLFWwindow* window, view_params* view) {
    char ext[1024];
    yc_split_path(view->imfilename.c_str(), 0, 0, ext);
    if (strcmp(ext, ".png")) {
        printf("supports only png screenshots");
        return;
    }

    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    ym_vector<unsigned char> pixels = ym_vector<unsigned char>(w * h * 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    ym_vector<unsigned char> line(w * 4);
    for (int j = 0; j < h / 2; j++) {
        memcpy(line.data(), pixels.data() + j * w * 4, w * 4);
        memcpy(pixels.data() + j * w * 4, pixels.data() + (h - 1 - j) * w * 4,
               w * 4);
        memcpy(pixels.data() + (h - 1 - j) * w * 4, line.data(), w * 4);
    }
    stbi_write_png(view->imfilename.c_str(), w, h, 4, pixels.data(), w * 4);
}

// text callback
void text_callback(GLFWwindow* window, unsigned int key) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    switch (key) {
        case '[': view->exposure -= 1; break;
        case ']': view->exposure += 1; break;
        case '{': view->gamma -= 0.1f; break;
        case '}': view->gamma += 0.1f; break;
        case 'w': view->wireframe = !view->wireframe; break;
        case 'e': view->edges = !view->edges; break;
        case 'b': view->background = (view->background + 1) % 4; break;
        case 's': save_screenshot(window, view); break;
        case 'c': view->camera_lights = !view->camera_lights; break;
        case 'C':
            view->cur_camera =
                (view->cur_camera + 1) % view->scene->cameras.size();
            break;
        case 't': {
            for (int i = 0; i < view->scene->shapes.size(); i++) {
                yo_shape* shape = &view->scene->shapes[i];
                ys_tesselate_stdshape(&shape->nelems, &shape->elem,
                                      shape->etype, &shape->nverts, &shape->pos,
                                      &shape->norm, &shape->texcoord,
                                      &shape->color, &shape->radius);
            }
        } break;
        default: printf("unsupported key\n"); break;
    }
}

void window_size_callback(GLFWwindow* window, int w, int h) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    view->w = w;
    view->h = h;
    yo_camera* cam = &view->scene->cameras[view->cur_camera];
    cam->width = cam->height * (float)view->w / (float)view->h;
}

void framebuffer_size_callback(GLFWwindow* window, int w, int h) {
    glViewport(0, 0, w, h);
}

void mouse_button_callback(GLFWwindow* window, int button, int action,
                           int mods) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    if (action == GLFW_RELEASE) {
        view->mouse_button = 0;
    } else if (button == GLFW_MOUSE_BUTTON_1 && !mods) {
        view->mouse_button = 1;
    } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_CONTROL)) {
        view->mouse_button = 2;
    } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_SHIFT)) {
        view->mouse_button = 3;
    } else if (button == GLFW_MOUSE_BUTTON_2) {
        view->mouse_button = 2;
    } else {
        view->mouse_button = 0;
    }
}

void mouse_pos_callback(GLFWwindow* window, double x, double y) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);

    ym_vec2f mouse_pos = {(float)x, (float)y};
    if (view->mouse_button) {
        float dolly = 0;
        ym_vec2f pan = ym_zero2f;
        ym_vec2f rotate = ym_zero2f;
        switch (view->mouse_button) {
            case 1: rotate = (mouse_pos - view->mouse_pos) / 100; break;
            case 2: dolly = (mouse_pos[0] - view->mouse_pos[0]) / 100.0f; break;
            case 3: pan = (mouse_pos - view->mouse_pos) / 100; break;
            default: break;
        }

        yo_camera* cam = &view->scene->cameras[view->cur_camera];
        ym_turntable(&cam->from, &cam->to, &cam->up, rotate, dolly, pan);
    }

    view->mouse_pos = mouse_pos;
}

void window_refresh_callback(GLFWwindow* window) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);

    char title[4096];
    sprintf(title, "mview | %s", view->filename.c_str());

    // begin frame
    ym_vec4f background = view->backgrounds[view->background];
    glClearColor(background[0], background[1], background[2], background[3]);
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // draw
    shade(view->scene, view->cur_camera, view->shade_prog,
          view->shade_txt.data(), view->exposure, view->gamma, view->wireframe,
          view->edges, view->camera_lights);

    // end frame
    glfwSwapBuffers(window);
}

// uiloop
void ui_loop(const ym_string& filename, const ym_string& imfilename,
             yo_scene* scene, int w, int h, int ns, bool camera_lights,
             float amb, bool no_ui) {
    // view data
    view_params* view =
        init_view_params(filename, imfilename, scene, w, h, camera_lights, amb);

    // glfw
    if (!glfwInit()) exit(EXIT_FAILURE);
    GLFWwindow* window = glfwCreateWindow(view->w, view->h, "yshade", 0, 0);
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, view);

    // callbacks
    glfwSetCharCallback(window, text_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_pos_callback);
    glfwSetWindowRefreshCallback(window, window_refresh_callback);

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) exit(EXIT_FAILURE);
#endif

    // init shade state
    view->shade_prog = yg_stdshader_make_program();
    view->shade_txt.resize(scene->textures.size());
    for (int i = 0; i < scene->textures.size(); i++) {
        yo_texture* txt = &scene->textures[i];
        view->shade_txt[i] =
            yg_make_texture(txt->pixels.data(), txt->width, txt->height,
                            txt->ncomp, false, true);
    }

    // ui loop
    while (!glfwWindowShouldClose(window)) {
        window_refresh_callback(window);

        if (no_ui) {
            printf("shading %s to %s\n", filename.c_str(), imfilename.c_str());
            save_screenshot(window, view);
            break;
        }

        glfwWaitEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    delete view;
}

// load scene and make intersect state
yo_scene* load_scene(const char* filename) {
    // load scenes and merge
    char ext[16];
    yc_split_path(filename, 0, 0, ext);
    yo_scene* scene = (strcmp(ext, ".objbin"))
                          ? yo_load_obj(filename, true, true)
                          : yo_load_objbin(filename, true);
    if (!scene) {
        printf("unable to load scene %s\n", filename);
        return 0;
    }

    // load textures
    yo_load_textures(scene, filename, 0);

    // check texture and return error if not found
    for (int t = 0; t < scene->textures.size(); t++) {
        yo_texture* txt = &scene->textures[t];
        if (txt->pixels.empty()) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->width = 1;
            txt->height = 1;
            txt->ncomp = 4;
            txt->pixels = {1, 1, 1, 1};
        }
    }

    // ensure normals
    for (int s = 0; s < scene->shapes.size(); s++) {
        yo_shape* shape = &scene->shapes[s];
        if (!shape->norm.empty()) continue;
        shape->norm.resize(shape->pos.size());
        ys_compute_normals(shape->nelems, shape->elem.data(), shape->etype,
                           shape->nverts, shape->pos.data(), shape->norm.data(),
                           true);
    }

    // make camera if not there
    if (!scene->cameras.size()) {
        // find scene bounds
        ym_range3f bbox = ym_invalid_range3f;
        for (int i = 0; i < scene->shapes.size(); i++) {
            yo_shape* shape = &scene->shapes[i];
            for (int j = 0; j < shape->nverts; j++) bbox += shape->pos[j];
        }
        ym_vec3f bbox_center = ym_rcenter(bbox);
        ym_vec3f bbox_size = ym_rsize(bbox);
        float bbox_msize = fmax(bbox_size.x, fmax(bbox_size.y, bbox_size.z));
        // create camera
        scene->cameras.resize(1);
        yo_camera* cam = &scene->cameras[0];
        // set up camera
        ym_vec3f camera_dir = {1, 0.4f, 1};
        cam->from = camera_dir * bbox_msize + bbox_center;
        cam->to = bbox_center;
        cam->up = ym_vec3f{0, 1, 0};
        cam->width = 16.0f / 9.0f;
        cam->height = 1;
        cam->aperture = 0;
    }

    return scene;
}

int main(int argc, const char** argv) {
    // command line
    yc_parser* parser = yc_init_parser(argc, argv, "view meshes");
    float amb = yc_parse_optf(parser, "--ambient", 0, "ambient factor", 0);
    bool camera_lights = yc_parse_optb(parser, "--camera_lights", "-c",
                                       "enable camera lights", false);
    int camera = yc_parse_opti(parser, "--camera", "-C", "camera", 0);
    bool no_ui = yc_parse_optb(parser, "--no-ui", 0, "runs offline", false);
    int samples = yc_parse_opti(parser, "--samples", "-s", "image samples", 1);
    float aspect =
        yc_parse_optf(parser, "--aspect", "-a", "image aspect", 16.0f / 9.0f);
    int res =
        yc_parse_opti(parser, "--resolution", "-r", "image resolution", 720);
    const char* imfilename =
        yc_parse_opts(parser, "--output", "-o", "image filename", "out.png");
    const char* filename =
        yc_parse_args(parser, "scene", "scene filename", 0, true);
    yc_done_parser(parser);

    // load scene
    yo_scene* scene = load_scene(filename);
    scene->cameras[camera].width = aspect * scene->cameras[camera].height;

    // start
    ui_loop(filename, imfilename, scene, (int)round(res * aspect), res, samples,
            camera_lights, amb, no_ui);

    // done
    yo_free_scene(scene);
    return EXIT_SUCCESS;
}
