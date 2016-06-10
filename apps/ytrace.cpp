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
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_glu.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#include "ext/ThreadPool1/ThreadPool.h"

#define MAX_BLOCKS_PER_STEP 4096

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

    // rendering params
    int ns = 0, nb = 0;
    int bs = 0, ss = 0;  // block size and step size
    yt_render_params params;
    ym_vector<ym_range2i> blocks;

    // current image step
    int cb = 0, cs = 0;

    // mouse
    int mouse_button = 0;
    ym_vec2f mouse_pos = ym_zero2f;

    // scene
    yo_scene* scene = nullptr;
    yt_scene* trace_scene = nullptr;
    int cur_camera = 0;

    // images
    ym_vector<ym_vec4f> img, buf, preview;

    // texture id
    int tid = 0;

    // scene sync
    bool scene_updated = false;

    // thread pool
    int nthreads = 0;
    ThreadPool* pool = nullptr;

    // stats
    bool stats_print_sps = false;
    double stats_ct = 0;
    double stats_sps = 0;
};

view_params* init_view_params(const ym_string& filename,
                              const ym_string& imfilename, yo_scene* scene,
                              yt_scene* trace_scene, int camera, int w, int h,
                              int ns, int bs, const yt_render_params& params,
                              int nthreads) {
    view_params* view = new view_params();

    view->filename = filename;
    view->imfilename = imfilename;

    view->w = w;
    view->h = h;

    view->ns = ns;
    view->bs = bs;
    view->params = params;
    view->ss = 8;

    view->scene = scene;
    view->trace_scene = trace_scene;

    view->cur_camera = camera;

    view->nthreads = nthreads;
    view->pool = new ThreadPool(nthreads);

    return view;
}

ym_vector<ym_range2i> make_image_blocks(int w, int h, int bs) {
    ym_vector<ym_range2i> blocks;
    for (int j = 0, cb = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs, cb++) {
            blocks.push_back({{i, j}, {ym_min(i + bs, w), ym_min(j + bs, h)}});
        }
    }
    return blocks;
}

void save_image(const ym_string& filename, ym_vec4f* pixels, int w, int h) {
    char ext[1024];
    yc_split_path(filename.c_str(), 0, 0, ext);
    if (!strcmp(ext, ".hdr")) {
        stbi_write_hdr(filename.c_str(), w, h, 4, (float*)pixels);
    } else if (!strcmp(ext, ".png")) {
        ym_vector<ym_vec4b> pixels_ub = ym_vector<ym_vec4b>(w * h);
        for (int i = 0; i < w * h; i++) {
            ym_vec4f p = pixels[i];
            ym_vec4b& pub = pixels_ub[i];
            pub[0] = ym_clamp((int)(256 * powf(p[0], 1 / 2.2f)), 0, 255);
            pub[1] = ym_clamp((int)(256 * powf(p[1], 1 / 2.2f)), 0, 255);
            pub[2] = ym_clamp((int)(256 * powf(p[2], 1 / 2.2f)), 0, 255);
            pub[3] = ym_clamp((int)(256 * p[3]), 0, 255);
        }
        stbi_write_png(filename.c_str(), w, h, 4, pixels_ub.data(), w * 4);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

// prepare image
void resize_images(view_params* view) {
    // allocate images
    view->img.resize(view->w * view->h);
    view->buf.resize(view->w * view->h);
    view->preview.resize(view->w / view->bs * view->h / view->bs);
    if (view->tid) yg_clear_texture(&view->tid);
    memset(view->img.data(), 0, sizeof(ym_vec4f) * view->w * view->h);

    // update texture
    // TODO: destroy preview textures
    view->tid = yg_make_texture((float*)view->img.data(), view->w, view->h, 4,
                                true, false);

    // create blocks
    view->blocks = make_image_blocks(view->w, view->h, view->bs);
    view->nb = (int)view->blocks.size();
}

// set renderer to restart
void render_restart(view_params* view) {
    view->cs = 0;
    view->cb = 0;
    view->stats_ct = 0;
}

// set lights
void update_trace_camera(view_params* view) {
    yo_camera* cam = &view->scene->cameras[view->cur_camera];
    ym_frame3f camera_frame = ym_lookat_xform3(cam->from, cam->to, cam->up);
    yt_set_camera(view->trace_scene, view->cur_camera, camera_frame, cam->width,
                  cam->height, cam->aperture, ym_dist(cam->from, cam->to));
}

// text callback
void text_callback(GLFWwindow* window, unsigned int key) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    switch (key) {
        case '[': view->exposure -= 1; break;
        case ']': view->exposure += 1; break;
        case '{': view->gamma -= 0.1f; break;
        case '}': view->gamma += 0.1f; break;
        case 's': {
            save_image(view->imfilename.c_str(), view->img.data(), view->w,
                       view->h);
        } break;
        case 'C': {
            view->cur_camera =
                (view->cur_camera + 1) % view->scene->cameras.size();
            update_trace_camera(view);
            view->scene_updated = true;
        } break;
        case 'b': view->background = (view->background + 1) % 4; break;
        case 'P': view->stats_print_sps = !view->stats_print_sps; break;
        default: printf("unsupported key\n"); break;
    }
}

void window_size_callback(GLFWwindow* window, int w, int h) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    view->w = w;
    view->h = h;
    yo_camera* cam = &view->scene->cameras[0];
    cam->width = cam->height * (float)view->w / (float)view->h;
    resize_images(view);
    view->scene_updated = true;
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

        view->scene_updated = true;
    }

    view->mouse_pos = mouse_pos;
}

void window_refresh_callback(GLFWwindow* window) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);

    char title[4096];
    sprintf(title, "ytrace | %dx%d | %d/%d  | %s", view->w, view->h, view->cs,
            view->ns, view->filename.c_str());
    glfwSetWindowTitle(window, title);

    ym_vec4f background = view->backgrounds[view->background];
    glClearColor(background[0], background[1], background[2], background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    yg_shade_image(view->tid, view->w, view->h, view->w, view->h, 0, 0, 1,
                   view->exposure, view->gamma);

    glfwSwapBuffers(window);
}

void render_preview(view_params* view) {
    yt_trace_image(view->trace_scene, view->cur_camera, view->preview.data(),
                   view->w / view->bs, view->h / view->bs, 1, view->params);
    for (int qj = 0; qj < view->h / view->bs; qj++) {
        for (int qi = 0; qi < view->w / view->bs; qi++) {
            const ym_vec4f& qpixel =
                view->preview[qj * (view->w / view->bs) + qi];
            for (int j = qj * view->bs;
                 j < ym_min((qj + 1) * view->bs, view->h); j++) {
                for (int i = qi * view->bs;
                     i < ym_min((qi + 1) * view->bs, view->w); i++) {
                    view->img[j * view->w + i] = qpixel;
                }
            }
        }
    }
    yg_update_texture(view->tid, (float*)view->img.data(), view->w, view->h, 4);
}

struct render_block_params {
    int cur_camera;
    ym_vec4f *img, *buf;
    int w, h, ns;
    ym_range2i window;
    ym_range1i samples;
    yt_render_params params;
    yt_scene* scene;
};

void accumulate_block(ym_vec4f* img, ym_vec4f* buf, int w, int h,
                      const ym_range2i& window, const ym_range1i& samples) {
    for (int jj = window.min.y; jj < window.max.y; jj++) {
        for (int ii = window.min.x; ii < window.max.x; ii++) {
            int off = (jj * w + ii);
            ym_vec4f *pixel = img + off, *bpixel = buf + off;
            if (samples.min) {
                *pixel = (*pixel * samples.min + *bpixel) / samples.max;
            } else {
                *pixel = *bpixel;
            }
        }
    }
}

void render_block_async(void* ctx) {
    render_block_params* params = (render_block_params*)ctx;
    yt_trace_block(params->scene, params->cur_camera,
                   (params->buf) ? (params->buf) : params->img, params->w,
                   params->h, params->ns, params->window, params->samples,
                   params->params);
    if (params->buf) {
        accumulate_block(params->img, params->buf, params->w, params->h,
                         params->window, params->samples);
    }
}

void render_step(view_params* view) {
    static render_block_params params[MAX_BLOCKS_PER_STEP];
    if (view->cs == view->ns) return;
    double timer_start = glfwGetTime();
    std::vector<std::future<void>> futures;
    for (int b = 0; view->cb < view->nb && b < view->ss; view->cb++, b++) {
        params[b].scene = view->trace_scene;
        params[b].img = view->img.data();
        params[b].buf = view->buf.data();
        params[b].w = view->w;
        params[b].h = view->h;
        params[b].ns = view->ns;
        params[b].window = view->blocks[view->cb];
        params[b].samples = {view->cs, view->cs + 1};
        params[b].params = view->params;
        futures.push_back(view->pool->enqueue(render_block_async, &params[b]));
    }
    for (int i = 0; i < futures.size(); i++) futures[i].wait();
    double timer_stop = glfwGetTime();
    view->stats_ct += timer_stop - timer_start;
    view->stats_sps = (double)(view->cs * view->w * view->h +
                               view->cb * view->bs * view->bs) /
                      view->stats_ct;
    yg_update_texture(view->tid, (float*)view->img.data(), view->w, view->h, 4);
    if (view->cb == view->nb) {
        view->cb = 0;
        view->cs++;
    }
}

// uiloop
void ui_loop(const ym_string& filename, const ym_string& imfilename,
             yo_scene* scene, yt_scene* trace_scene, int camera, int w, int h,
             int ns, int bs, bool camera_lights, const yt_render_params& params,
             int nthreads) {
    // view data
    view_params* view =
        init_view_params(filename, imfilename, scene, trace_scene, camera, w, h,
                         ns, bs, params, nthreads);

    // glfw
    if (!glfwInit()) exit(EXIT_FAILURE);
    GLFWwindow* window = glfwCreateWindow(view->w, view->h, "mtrace", 0, 0);
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, view);

    // callbacks
    glfwSetCharCallback(window, text_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_pos_callback);
    glfwSetWindowRefreshCallback(window, window_refresh_callback);

#ifndef __APPLE__
    // init gl extensions
    if (glewInit() != GLEW_OK) exit(EXIT_FAILURE);
#endif

    // prepare images
    resize_images(view);
    view->scene_updated = true;

    // ui loop
    while (!glfwWindowShouldClose(window)) {
        if (view->scene_updated) {
            update_trace_camera(view);
            render_preview(view);
            render_restart(view);
            view->scene_updated = false;
        } else {
            render_step(view);
            if (view->stats_print_sps) printf("%.0lf\n", view->stats_sps);
        }

        window_refresh_callback(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    delete view->pool;

    delete view;
}

bool intersect_first(const void* ctx, const ym_ray3f& ray, float* ray_t,
                     int* sid, int* eid, ym_vec2f* euv) {
    const yb_scene* scene_bvh = (const yb_scene*)ctx;
    return yb_intersect_first(scene_bvh, ray, ray_t, sid, eid, euv);
}

bool intersect_any(const void* ctx, const ym_ray3f& ray) {
    const yb_scene* scene_bvh = (const yb_scene*)ctx;
    return yb_intersect_any(scene_bvh, ray);
}

yb_scene* make_bvh(const yo_scene* scene) {
    yb_scene* scene_bvh = yb_init_scene((int)scene->shapes.size());
    for (int sid = 0; sid < scene->shapes.size(); sid++) {
        const yo_shape* shape = &scene->shapes[sid];
        yb_set_shape(scene_bvh, sid, shape->xform, shape->nelems,
                     shape->elem.data(), shape->etype, shape->nverts,
                     shape->pos.data(), shape->radius.data());
    }
    yb_build_bvh(scene_bvh, 0);
    return scene_bvh;
}

yt_scene* init_trace_cb(const yo_scene* scene, yb_scene* scene_bvh,
                        int camera) {
    yt_scene* trace_scene =
        yt_make_scene((int)scene->cameras.size(), (int)scene->shapes.size(),
                      (int)scene->materials.size(), (int)scene->textures.size(),
                      (int)scene->envs.size());

    yt_set_intersection(trace_scene, scene_bvh, intersect_first, intersect_any);

    for (int sid = 0; sid < scene->envs.size(); sid++) {
        const yo_camera* cam = &scene->cameras[sid];
        ym_frame3f camera_frame = ym_lookat_xform3(cam->from, cam->to, cam->up);
        yt_set_camera(trace_scene, sid, camera_frame, cam->width, cam->height,
                      cam->aperture, ym_dist(cam->from, cam->to));
    }

    for (int sid = 0; sid < scene->envs.size(); sid++) {
        const yo_env* env = &scene->envs[sid];
        ym_frame3f env_frame = ym_lookat_xform3(env->from, env->to, env->up);
        const yo_material* mat = &scene->materials[env->matid];
        yt_set_env(trace_scene, sid, env_frame, (ym_vec3f)mat->ke,
                   mat->ke_txtid);
    }

    for (int sid = 0; sid < scene->shapes.size(); sid++) {
        const yo_shape* shape = &scene->shapes[sid];
        ym_frame3f shape_frame = (ym_frame3f)(ym_mat4f)shape->xform;
        yt_set_shape(trace_scene, sid, shape->matid, shape_frame, shape->nelems,
                     shape->elem.data(), shape->etype, shape->nverts,
                     shape->pos.data(), shape->norm.data(),
                     shape->texcoord.data(), 0, 0);
    }

    for (int mid = 0; mid < scene->materials.size(); mid++) {
        const yo_material* mat = &scene->materials[mid];
        yt_set_material(trace_scene, mid, (ym_vec3f)mat->ke, (ym_vec3f)mat->kd,
                        (ym_vec3f)mat->ks,
                        yt_specular_exponent_to_roughness(mat->ns), ym_zero3f,
                        ym_zero3f, mat->ke_txtid, mat->kd_txtid, mat->ks_txtid,
                        mat->ns_txtid, false);
    }

    for (int tid = 0; tid < scene->textures.size(); tid++) {
        const yo_texture* txt = &scene->textures[tid];
        yt_set_texture(trace_scene, tid, txt->pixels.data(), txt->width,
                       txt->height, txt->ncomp);
    }

    yt_init_lights(trace_scene);

    return trace_scene;
}

void render(yt_scene* trace_scene, const char* filename, const char* imfilename,
            int cid, int w, int h, int ns, const yt_render_params& params,
            int bs, int nthreads) {
    ym_vector<ym_vec4f> pixels = ym_vector<ym_vec4f>(w * h);
    ym_vector<ym_vec4f> buf = ym_vector<ym_vec4f>(w * h);
    ym_vector<ym_range2i> blocks = make_image_blocks(w, h, bs);
    ym_vector<render_block_params> block_params =
        ym_vector<render_block_params>(blocks.size());
    ThreadPool pool(nthreads);
    std::vector<std::future<void>> futures;
    printf("tracing %s to %s\n", filename, imfilename);
    printf("rendering ...");
    fflush(stdout);
    for (int s = 0; s < ns; s++) {
        printf("\rrendering sample %d/%d", s + 1, ns);
        fflush(stdout);
        for (int i = 0; i < blocks.size(); i++) {
            block_params[i] = render_block_params{
                cid, pixels.data(), buf.data(), w,      h,
                ns,  blocks[i],     {s, s + 1}, params, trace_scene};
            futures.emplace_back(
                pool.enqueue(render_block_async, &block_params[i]));
        }
        for (int i = 0; i < futures.size(); i++) futures[i].wait();
    }
    printf("\rrendering done\n");
    fflush(stdout);
    save_image(imfilename, pixels.data(), w, h);
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

    // ensure radius is necessary
    for (int s = 0; s < scene->shapes.size(); s++) {
        yo_shape* shape = &scene->shapes[s];
        if (shape->etype != yo_etype_point && shape->etype != yo_etype_line)
            continue;
        if (!shape->radius.empty()) continue;
        shape->radius.resize(shape->pos.size());
        for (int i = 0; i < shape->nverts; i++) shape->radius[i] = 0.001f;
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
    const char* rtype_names[] = {"default", "uniform", "stratified", "cmjs", 0};
    const char* stype_names[] = {"default", "eye", "direct", "path", 0};
    // params
    yc_parser* parser = yc_init_parser(argc, argv, "trace meshes");
    int rtype =
        yc_parse_opte(parser, "--random", 0, "random type", 0, rtype_names);
    int stype = yc_parse_opte(parser, "--integrator", "-i", "integrator type",
                              0, stype_names);
    float amb = yc_parse_optf(parser, "--ambient", 0, "ambient factor", 0);
    bool camera_lights = yc_parse_optb(parser, "--camera_lights", "-c",
                                       "enable camera lights", false);
    int camera = yc_parse_opti(parser, "--camera", "-C", "camera", 0);
    bool no_ui = yc_parse_optb(parser, "--no-ui", 0, "runs offline", false);
    int nthreads = yc_parse_opti(parser, "--threads", "-t",
                                 "number of threads [0 for default]", 0);
    int block_size = yc_parse_opti(parser, "--block_size", 0, "block size", 32);
    int samples =
        yc_parse_opti(parser, "--samples", "-s", "image samples", 256);
    float aspect =
        yc_parse_optf(parser, "--aspect", "-a", "image aspect", 16.0f / 9.0f);
    int res =
        yc_parse_opti(parser, "--resolution", "-r", "image resolution", 720);
    const char* imfilename =
        yc_parse_opts(parser, "--output", "-o", "image filename", "out.hdr");
    const char* filename =
        yc_parse_args(parser, "scene", "scene filename", 0, true);
    yc_done_parser(parser);

    // setting up multithreading
    if (!nthreads) nthreads = std::thread::hardware_concurrency();

    // loading scene
    yo_scene* scene = load_scene(filename);
    scene->cameras[camera].width = aspect * scene->cameras[camera].height;

    // preparing raytracer
    yb_scene* scene_bvh = make_bvh(scene);
    yt_scene* trace_scene = init_trace_cb(scene, scene_bvh, camera);
    yt_render_params params;
    params.stype = (camera_lights) ? yt_stype_eyelight : stype;
    params.rtype = rtype;
    params.amb = {amb, amb, amb};

    // launching renderer
    if (no_ui) {
        render(trace_scene, filename, imfilename, camera,
               (int)round(res * aspect), res, samples, params, block_size,
               nthreads);
    } else {
        ui_loop(filename, imfilename, scene, trace_scene, camera,
                (int)round(res * aspect), res, samples, block_size,
                camera_lights, params, nthreads);
    }

    // done
    // TODO: free bvh
    yt_free_scene(trace_scene);
    yo_free_scene(scene);
    return EXIT_SUCCESS;
}
