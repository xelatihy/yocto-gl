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

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"

#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

// view data
struct view_params {
    // filenames
    const char* filename;
    const char* imfilename;

    // window data
    int w, h;
    float exposure, gamma;
    float* background;
    float backgrounds[16];

    // drawing
    bool wireframe;
    bool edges;

    // mouse state
    int mouse_button;
    double mouse_x, mouse_y;

    // scene
    yo_scene* scene;
    float amb[3];
    int cur_camera;

    // shading
    int shade_prog;
    int shade_ntxt;
    int* shade_txt;

    // lighting
    bool camera_lights;
};
typedef struct view_params view_params;

view_params*
init_view_params(const char* filename, const char* imfilename, yo_scene* scene,
                 int w, int h, bool camera_lights, float amb) {
    view_params* view = (view_params*)calloc(1, sizeof(view_params));

    view->filename = filename;
    view->imfilename = imfilename;

    view->w = w;
    view->h = h;
    view->exposure = 0;
    view->gamma = 2.2f;
    view->background = view->backgrounds;

    float backgrounds[16] = { 0,   0,   0,   0, 0.18, 0.18, 0.18, 0,
                              0.5, 0.5, 0.5, 0, 1,    1,    1,    0 };
    memcpy(view->backgrounds, backgrounds, sizeof(view->backgrounds));

    view->wireframe = false;
    view->edges = false;

    view->scene = scene;
    view->amb[0] = view->amb[1] = view->amb[2] = amb;

    view->camera_lights = camera_lights;

    return view;
}

void
turntable(ym_vec3f* from, ym_vec3f* to, ym_vec3f* up, const float rotate[2],
          float dolly, const float pan[2]) {
    // rotate if necessary
    if (rotate && (rotate[0] || rotate[1])) {
        ym_vec3f z = ym_normalize3f(ym_sub3f(*to, *from));
        float lz = ym_dist3f(*to, *from);
        float phi = atan2f(z.z, z.x) + rotate[0];
        float theta = acosf(z.y) + rotate[1];
        theta = fmax(0.001f, fmin(theta, ym_pif - 0.001f));
        ym_vec3f nz = { sinf(theta) * cosf(phi) * lz, cosf(theta) * lz,
                        sinf(theta) * sinf(phi) * lz };
        *from = ym_sub3f(*to, nz);
    }

    // dolly if necessary
    if (dolly) {
        ym_vec3f z = ym_normalize3f(ym_sub3f(*to, *from));
        float lz = fmaxf(0.001f, ym_dist3f(*to, *from) * (1 + dolly));
        z = ym_smul3f(z, lz);
        *from = ym_sub3f(*to, z);
    }

    // pan if necessary
    if (pan) {
        ym_vec3f z = ym_normalize3f(ym_sub3f(*to, *from));
        ym_vec3f x = ym_normalize3f(ym_cross3f(*up, z));
        ym_vec3f y = ym_normalize3f(ym_cross3f(z, x));
        ym_vec3f t = { pan[0] * x.x + pan[1] * y.x, pan[0] * x.y + pan[1] * y.y,
                       pan[0] * x.z + pan[1] * y.z };
        *from = ym_sum3f(*from, t);
        *to = ym_sum3f(*to, t);
    }
}

void
shade(yo_scene* scene, int cur_camera, int prog, int* txt, float exposure,
      float gamma_, bool wireframe, bool edges, bool camera_lights) {
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    yo_camera* cam = &scene->cameras[cur_camera];
    ym_mat4f camera_xform = ym_lookat_xform4f(
        *(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to, *(ym_vec3f*)cam->up);
    ym_mat4f camera_view = ym_lookat_view4f(
        *(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to, *(ym_vec3f*)cam->up);
    ym_mat4f camera_proj = ym_perspective4f(
        2 * atanf(cam->height / 2), cam->width / cam->height, 0.1f, 10000);

    yg_stdshader_begin_frame(prog, camera_lights, exposure, gamma_,
                             camera_xform.m, camera_view.m, camera_proj.m);

    if (!camera_lights) {
        int nlights = 0;
        ym_vec3f light_pos[16], light_ke[16];
        int light_type[16];
        for (int i = 0; i < scene->nshapes && nlights < 16; i++) {
            yo_shape* shape = &scene->shapes[i];
            if (shape->etype != yo_etype_point) continue;
            if (shape->matid < 0) continue;
            yo_material* mat = &scene->materials[shape->matid];
            if (!mat->ke[0] && !mat->ke[1] && !mat->ke[2]) continue;
            for (int j = 0; j < shape->nverts && nlights < 16; j++) {
                light_pos[nlights] = *(ym_vec3f*)(shape->pos + j * 3);
                if (shape->xform) {
                    light_pos[nlights] = ym_transform_point3f(
                        *(ym_mat4f*)shape->xform, light_pos[nlights]);
                }
                light_ke[nlights] = *(ym_vec3f*)mat->ke;
                light_type[nlights] = 0;
                nlights++;
            }
        }
        ym_vec3f zero = { 0, 0, 0 };
        yg_stdshader_set_lights(prog, &zero.x, nlights, (float*)light_pos,
                                (float*)light_ke, light_type);
    }

    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = &scene->shapes[i];

        ym_mat4f shape_xform =
            (shape->xform) ? *(ym_mat4f*)shape->xform : ym_identity4f();
        yg_stdshader_begin_shape(prog, shape_xform.m);

        if (shape->matid >= 0) {
            yo_material* mat = &scene->materials[shape->matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yg_stdshader_set_material(
                prog, shape->etype, mat->ke, mat->kd, mat->ks,
                yg_specular_exponent_to_roughness(mat->ns),
                __txt(mat->ke_txtid), __txt(mat->kd_txtid),
                __txt(mat->ks_txtid), __txt(mat->ns_txtid), false);
        } else {
            float kd[3] = { 0.8, 0.8, 0.8 }, zero[3] = { 0, 0, 0 };
            yg_stdshader_set_material(prog, shape->etype, zero, kd, zero, 0, 0,
                                      0, 0, 0, false);
        }

        yg_stdshader_set_vert(prog, shape->pos, shape->norm, shape->texcoord,
                              0);

        yg_stdshader_draw_elem(prog, shape->nelems, shape->elem, shape->etype);

        if (edges && !wireframe) {
            float zero[3] = { 0, 0, 0 };
            yg_stdshader_set_material(prog, shape->etype, zero, zero, zero, 0,
                                      0, 0, 0, 0, false);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yg_stdshader_draw_elem(prog, shape->nelems, shape->elem,
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

void
save_screenshot(GLFWwindow* window, view_params* view) {
    char ext[1024];
    yc_split_path(view->imfilename, 0, 0, ext);
    if (strcmp(ext, ".png")) {
        printf("supports only png screenshots");
        return;
    }

    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    unsigned char* pixels =
        (unsigned char*)calloc(w * h * 4, sizeof(unsigned char));
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    unsigned char* line = (unsigned char*)calloc(w * 4, sizeof(unsigned char));
    for (int j = 0; j < h / 2; j++) {
        memcpy(line, pixels + j * w * 4, w * 4);
        memcpy(pixels + j * w * 4, pixels + (h - 1 - j) * w * 4, w * 4);
        memcpy(pixels + (h - 1 - j) * w * 4, line, w * 4);
    }
    stbi_write_png(view->imfilename, w, h, 4, pixels, w * 4);
    free(line);
    free(pixels);
}

// text callback
void
text_callback(GLFWwindow* window, unsigned int key) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    switch (key) {
        case '[': view->exposure -= 1; break;
        case ']': view->exposure += 1; break;
        case '{': view->gamma -= 0.1; break;
        case '}': view->gamma += 0.1; break;
        case 'w': view->wireframe = !view->wireframe; break;
        case 'e': view->edges = !view->edges; break;
        case 'b':
            view->background += 4;
            if (view->background - view->backgrounds >= 16)
                view->background = view->backgrounds;
            break;
        case 's': save_screenshot(window, view); break;
        case 'c': view->camera_lights = !view->camera_lights; break;
        case 'C':
            view->cur_camera = (view->cur_camera + 1) % view->scene->ncameras;
            break;
        case 't': {
            for (int i = 0; i < view->scene->nshapes; i++) {
                yo_shape* shape = &view->scene->shapes[i];
                float** vert[5] = { &shape->pos, &shape->norm, &shape->texcoord,
                                    &shape->color, &shape->radius };
                int vsize[5] = { 3, 3, 2, 3, 1 };
                if (shape->etype == yo_etype_triangle ||
                    shape->etype == yo_etype_quad ||
                    shape->etype == yo_etype_line) {
                    ys_tesselate_shape(shape->etype, 5, vsize, &shape->nelems,
                                       &shape->elem, &shape->nverts, vert);
                }
            }
        } break;
        default: printf("unsupported key\n"); break;
    }
}

void
window_size_callback(GLFWwindow* window, int w, int h) {
    view_params* view = glfwGetWindowUserPointer(window);
    view->w = w;
    view->h = h;
    yo_camera* cam = &view->scene->cameras[view->cur_camera];
    cam->width = cam->height * (float)view->w / (float)view->h;
}

void
framebuffer_size_callback(GLFWwindow* window, int w, int h) {
    glViewport(0, 0, w, h);
}

void
mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    view_params* view = glfwGetWindowUserPointer(window);
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

void
mouse_pos_callback(GLFWwindow* window, double x, double y) {
    view_params* view = glfwGetWindowUserPointer(window);

    if (view->mouse_button) {
        float dolly = 0, pan[2] = { 0, 0 }, rotate[2] = { 0, 0 };
        switch (view->mouse_button) {
            case 1:
                rotate[0] = (x - view->mouse_x) / 100.0f;
                rotate[1] = (y - view->mouse_y) / 100.0f;
                break;
            case 2: dolly = (x - view->mouse_x) / 100.0f; break;
            case 3:
                pan[0] = (x - view->mouse_x) / 100.0f;
                pan[1] = (y - view->mouse_y) / 100.0f;
                break;
            default: break;
        }

        yo_camera* cam = &view->scene->cameras[view->cur_camera];
        turntable((ym_vec3f*)cam->from, (ym_vec3f*)cam->to, (ym_vec3f*)cam->up,
                  rotate, dolly, pan);
    }

    view->mouse_x = x;
    view->mouse_y = y;
}

void
window_refresh_callback(GLFWwindow* window) {
    view_params* view = glfwGetWindowUserPointer(window);

    char title[4096];
    sprintf(title, "mview | %s", view->filename);

    // begin frame
    glClearColor(view->background[0], view->background[1], view->background[2],
                 view->background[3]);
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // draw
    shade(view->scene, view->cur_camera, view->shade_prog, view->shade_txt,
          view->exposure, view->gamma, view->wireframe, view->edges,
          view->camera_lights);

    // end frame
    glfwSwapBuffers(window);
}

// uiloop
void
ui_loop(const char* filename, const char* imfilename, yo_scene* scene, int w,
        int h, int ns, bool camera_lights, float amb, bool no_ui) {
    // view data
    view_params* view =
        init_view_params(filename, imfilename, scene, w, h, camera_lights, amb);

    // glfw
    if (!glfwInit()) exit(EXIT_FAILURE);
    GLFWwindow* window = glfwCreateWindow(view->w, view->h, "mview", 0, 0);
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, view);

    // callbacks
    glfwSetCharCallback(window, text_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_pos_callback);
    glfwSetWindowRefreshCallback(window, window_refresh_callback);

    // initialize glew
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) exit(EXIT_FAILURE);

    // init shade state
    view->shade_prog = yg_stdshader_make_program();
    view->shade_ntxt = scene->ntextures;
    view->shade_txt = (int*)calloc(scene->ntextures, sizeof(int));
    for (int i = 0; i < scene->ntextures; i++) {
        yo_texture* txt = &scene->textures[i];
        view->shade_txt[i] = yg_make_texture(
            txt->pixels, txt->width, txt->height, txt->ncomp, false, true);
    }

    // ui loop
    while (!glfwWindowShouldClose(window)) {
        window_refresh_callback(window);

        if (no_ui) {
            printf("shading %s to %s\n", filename, imfilename);
            save_screenshot(window, view);
            break;
        }

        glfwWaitEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    free(view);
}

// load scene and make intersect state
yo_scene*
load_scene(const char* filename, bool triangulate) {
    // load scenes and merge
    char ext[16];
    yc_split_path(filename, 0, 0, ext);
    yo_scene* scene = (strcmp(ext, ".objbin"))
                          ? yo_load_obj(filename, triangulate, true)
                          : yo_load_objbin(filename, true);
    if (!scene) return 0;

    // load textures
    yo_load_textures(scene, filename, 0);

    // ensure normals
    for (int s = 0; s < scene->nshapes; s++) {
        yo_shape* shape = &scene->shapes[s];
        if (shape->norm) continue;
        shape->norm =
            ys_compute_normals(shape->nelems, shape->elem, shape->etype,
                               shape->nverts, shape->pos, shape->norm, true);
    }

    // make camera if not there
    if (!scene->ncameras) {
        // find scene bounds
        ym_range3f bbox = ym_rinit3f();
        for (int i = 0; i < scene->nshapes; i++) {
            yo_shape* shape = &scene->shapes[i];
            for (int j = 0; j < shape->nverts; j++) {
                bbox = ym_rexpand3f(bbox, ((ym_vec3f*)shape->pos)[j]);
            }
        }
        ym_vec3f bbox_center = ym_rcenter3f(bbox);
        ym_vec3f bbox_size = ym_rsize3f(bbox);
        float bbox_msize = fmax(bbox_size.x, fmax(bbox_size.y, bbox_size.z));
        // create camera
        scene->ncameras = 1;
        scene->cameras = (yo_camera*)calloc(1, sizeof(yo_camera));
        yo_camera* cam = &scene->cameras[0];
        // set up camera
        ym_vec3f camera_dir = { 1, 0.4f, 1 };
        *(ym_vec3f*)cam->from = ym_smul3f(camera_dir, bbox_msize);
        *(ym_vec3f*)cam->from = ym_sum3f(*(ym_vec3f*)cam->from, bbox_center);
        *(ym_vec3f*)cam->to = bbox_center;
        *(ym_vec3f*)cam->up = (ym_vec3f){ 0, 1, 0 };
        cam->width = 16.0f / 9.0f;
        cam->height = 1;
        cam->aperture = 0;
    }

    return scene;
}

int
main(int argc, const char** argv) {
    // command line
    yc_parser* parser = yc_init_parser(argc, argv, "view meshes");
    float amb = yc_parse_optf(parser, "--ambient", 0, "ambient factor", 0);
    bool camera_lights = yc_parse_optb(parser, "--camera_lights", "-c",
                                       "enable camera lights", false);
    int camera = yc_parse_opti(parser, "--camera", "-C", "camera", 0);
    bool no_ui = yc_parse_optb(parser, "--no-ui", 0, "runs offline", false);
    bool triangulate =
        yc_parse_optb(parser, "--triangulate", 0, "triangulate input", false);
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
    yo_scene* scene = load_scene(filename, triangulate);
    scene->cameras[camera].width = aspect * scene->cameras[camera].height;

    // start
    ui_loop(filename, imfilename, scene, round(res * aspect), res, samples,
            camera_lights, amb, no_ui);

    // done
    yo_free_scene(scene);
    return EXIT_SUCCESS;
}
