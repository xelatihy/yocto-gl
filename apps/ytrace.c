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
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"

#ifndef YA_NOGL
// clang-format off
#include "../yocto/yocto_glu.h"
#include <GLFW/glfw3.h>
// clang-format on
#endif

#ifndef _MSC_VER
#include "ext/thpool.c"
#include "ext/thpool.h"
#else
typedef int threadpool;
static inline threadpool
thpool_init(int n) {
    return 0;
}
static inline void
thpool_destroy(threadpool pool) {}
static inline void
thpool_wait(threadpool pool) {}
static inline void
thpool_add_work(threadpool pool, void* (*fn)(void*), void* ctx) {
    fn(ctx);
}

#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#ifndef YA_NOGL

#include <sys/time.h>

#define MAX_BLOCKS_PER_STEP 4096

// view data
struct view_params {
    // filenames
    const char* filename;
    const char* imfilename;

    // window data
    int w, h, nc;
    float exposure, gamma;
    float* background;
    float backgrounds[16];

    // rendering params
    int ns, nb;
    int bs, ss;  // block size and step size
    int stype, rtype;
    float amb[3];
    int* blocks;

    // current image step
    int cb, cs;

    // mouse
    int mouse_button;
    double mouse_x, mouse_y;

    // scene
    yo_scene* scene;
    yt_scene* trace_scene;
    int cur_camera;

    // lighting
    bool camera_lights;

    // images
    float *img, *buf, *preview;

    // texture id
    int tid;

    // scene sync
    bool scene_updated;

    // thread pool
    int nthreads;
    threadpool pool;

    // stats
    bool stats_print_sps;
    uint64_t stats_ct;
    double stats_sps;
};
typedef struct view_params view_params;

view_params*
init_view_params(const char* filename, const char* imfilename, yo_scene* scene,
                 yt_scene* trace_scene, int w, int h, int ns, int bs,
                 bool camera_lights, int stype, int rtype, float amb) {
    view_params* view = (view_params*)calloc(1, sizeof(view_params));

    view->filename = filename;
    view->imfilename = imfilename;

    view->w = w;
    view->h = h;
    view->nc = 4;
    view->exposure = 0;
    view->gamma = 2.2f;
    view->background = view->backgrounds;

    float backgrounds[16] = { 0,   0,   0,   0, 0.18, 0.18, 0.18, 0,
                              0.5, 0.5, 0.5, 0, 1,    1,    1,    0 };
    memcpy(view->backgrounds, backgrounds, sizeof(view->backgrounds));

    view->ns = ns;
    view->bs = bs;
    view->stype = stype;
    view->rtype = rtype;
    view->amb[0] = view->amb[1] = view->amb[2] = amb;
    view->ss = 8;

    view->scene = scene;
    view->trace_scene = trace_scene;

    view->camera_lights = camera_lights;

    view->nthreads = 4;
    view->pool = thpool_init(view->nthreads);

    return view;
}

#endif

int*
make_image_blocks(int w, int h, int bs, int* nb) {
    *nb = 0;
    for (int j = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs) (*nb)++;
    }
    int* blocks = (int*)calloc(*nb, 4 * sizeof(int));
    for (int j = 0, cb = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs, cb++) {
            blocks[cb * 4 + 0] = i;
            blocks[cb * 4 + 1] = ym_min(i + bs, w);
            blocks[cb * 4 + 2] = j;
            blocks[cb * 4 + 3] = ym_min(j + bs, h);
        }
    }
    return blocks;
}

void
save_image(const char* filename, float* pixels, int w, int h, int nc) {
    char ext[1024];
    yc_split_path(filename, 0, 0, ext);
    if (!strcmp(ext, ".hdr")) {
        stbi_write_hdr(filename, w, h, nc, pixels);
    } else if (!strcmp(ext, ".png")) {
        unsigned char* pixels_ub = (unsigned char*)malloc(w * h * nc);
        for (int i = 0; i < w * h; i++) {
            float* p = pixels + i * nc;
            unsigned char* pub = pixels_ub + i * nc;
            pub[0] = ym_clamp(256 * powf(p[0], 1 / 2.2f), 0, 255);
            pub[1] = ym_clamp(256 * powf(p[1], 1 / 2.2f), 0, 255);
            pub[2] = ym_clamp(256 * powf(p[2], 1 / 2.2f), 0, 255);
            if (nc == 4) pub[3] = ym_clamp(256 * p[3], 0, 255);
        }
        stbi_write_png(filename, w, h, nc, pixels_ub, w * nc);
        free(pixels_ub);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

#ifndef YA_NOGL

// prepare image
void
resize_images(view_params* view) {
    // allocate images
    view->img = (float*)realloc(view->img,
                                sizeof(float) * view->nc * view->w * view->h);
    view->buf = (float*)realloc(view->buf,
                                sizeof(float) * view->nc * view->w * view->h);
    view->preview =
        (float*)realloc(view->preview, sizeof(float) * view->nc * view->w /
                                           view->bs * view->h / view->bs);
    if (view->tid) yg_clear_texture(&view->tid);
    memset(view->img, 0, sizeof(float) * view->w * view->h * view->nc);

    // update texture
    // TODO: destroy preview textures
    view->tid =
        yg_make_texture(view->img, view->w, view->h, view->nc, true, false);

    // create blocks
    if (view->blocks) free(view->blocks);
    view->blocks = make_image_blocks(view->w, view->h, view->bs, &view->nb);
}

// set renderer to restart
void
render_restart(view_params* view) {
    view->cs = 0;
    view->cb = 0;
    view->stats_ct = 0;
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

// set lights
void
update_trace_camera(view_params* view) {
    yo_camera* cam = &view->scene->cameras[view->cur_camera];
    ym_mat4f camera_xform = ym_lookat_xform4f(
        *(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to, *(ym_vec3f*)cam->up);
    yt_set_camera(view->trace_scene, camera_xform.m, cam->width, cam->height,
                  cam->aperture,
                  ym_dist3f(*(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to));
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
        case 's': {
            save_image(view->imfilename, view->img, view->w, view->h, view->nc);
        } break;
        case 'c': {
            view->camera_lights = !view->camera_lights;
            if (view->camera_lights) {
                yt_set_rendering_params(view->trace_scene, yt_stype_eyelight,
                                        view->rtype, view->amb);
            } else {
                yt_set_rendering_params(view->trace_scene, view->stype,
                                        view->rtype, view->amb);
            }
            view->scene_updated = true;
        } break;
        case 'C': {
            view->cur_camera = (view->cur_camera + 1) % view->scene->ncameras;
            update_trace_camera(view);
            view->scene_updated = true;
        } break;
        case 'b':
            view->background += 4;
            if (view->background - view->backgrounds >= 16)
                view->background = view->backgrounds;
            break;
        case 'P': view->stats_print_sps = !view->stats_print_sps; break;
        default: printf("unsupported key\n"); break;
    }
}

void
window_size_callback(GLFWwindow* window, int w, int h) {
    view_params* view = glfwGetWindowUserPointer(window);
    view->w = w;
    view->h = h;
    yo_camera* cam = &view->scene->cameras[0];
    cam->width = cam->height * (float)view->w / (float)view->h;
    resize_images(view);
    view->scene_updated = true;
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
        view->scene_updated = true;
    }

    view->mouse_x = x;
    view->mouse_y = y;
}

void
window_refresh_callback(GLFWwindow* window) {
    view_params* view = glfwGetWindowUserPointer(window);

    char title[4096];
    sprintf(title, "mtrace | %dx%d | %d/%d  | %s", view->w, view->h, view->cs,
            view->ns, view->filename);
    glfwSetWindowTitle(window, title);

    glClearColor(view->background[0], view->background[1], view->background[2],
                 view->background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    yg_shade_image(view->tid, view->w, view->h, view->w, view->h, 0, 0, 1,
                   view->exposure, view->gamma);

    glfwSwapBuffers(window);
}

void
render_preview(view_params* view) {
    for (int qj = 0; qj < view->h / view->bs; qj++) {
        int block[6] = { 0, view->w / view->bs, qj, qj + 1, 0, 1 };
        yt_trace_block(view->trace_scene, view->preview, view->w / view->bs,
                       view->h / view->bs, 1, block);
        for (int qi = 0; qi < view->w / view->bs; qi++) {
            float* qpixel =
                view->preview + (qj * (view->w / view->bs) + qi) * 4;
            for (int j = qj * view->bs;
                 j < ym_min((qj + 1) * view->bs, view->h); j++) {
                for (int i = qi * view->bs;
                     i < ym_min((qi + 1) * view->bs, view->w); i++) {
                    float* pixel = view->img + (j * view->w + i) * 4;
                    for (int i = 0; i < 4; i++) pixel[i] = qpixel[i];
                }
            }
        }
    }
    yg_update_texture(view->tid, view->img, view->w, view->h, view->nc);
}

#endif

struct render_block_params {
    float *img, *buf;
    int w, h, ns;
    int block[6];
    yt_scene* scene;
};
typedef struct render_block_params render_block_params;

void
accumulate_block(float* img, float* buf, int w, int h, int block[6]) {
    for (int jj = block[2]; jj < block[3]; jj++) {
        for (int ii = block[0]; ii < block[1]; ii++) {
            int off = (jj * w + ii) * 4;
            float *pixel = img + off, *bpixel = buf + off;
            if (block[4]) {
                for (int i = 0; i < 4; i++)
                    pixel[i] = (pixel[i] * block[4] + bpixel[i]) / block[5];
            } else {
                for (int i = 0; i < 4; i++) pixel[i] = bpixel[i];
            }
        }
    }
}

void*
render_block_async(void* ctx) {
    render_block_params* params = (render_block_params*)ctx;
    yt_trace_block(params->scene, (params->buf) ? (params->buf) : params->img,
                   params->w, params->h, params->ns, params->block);
    if (params->buf)
        accumulate_block(params->img, params->buf, params->w, params->h,
                         params->block);
    return 0;
}

#ifndef YA_NOGL

uint64_t
get_timeus() {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (uint64_t)tv.tv_sec * 1000000ull + (uint64_t)tv.tv_usec;
}

void
render_step(view_params* view) {
    static render_block_params params[MAX_BLOCKS_PER_STEP];
    if (view->cs == view->ns) return;
    uint64_t timer_start = get_timeus();
    for (int b = 0; view->cb < view->nb && b < view->ss; view->cb++, b++) {
        int s = view->cs;
        int* imb = view->blocks + view->cb * 4;
        int block[6] = { imb[0], imb[1], imb[2], imb[3], s, s + 1 };
        if (view->nthreads) {
            params[b].scene = view->trace_scene;
            params[b].img = view->img;
            params[b].buf = view->buf;
            params[b].w = view->w;
            params[b].h = view->h;
            params[b].ns = view->ns;
            memcpy(params[b].block, block, sizeof(int) * 6);
            thpool_add_work(view->pool, render_block_async, &params[b]);
        } else {
            yt_trace_block(view->trace_scene, view->buf, view->w, view->h,
                           view->ns, block);
            accumulate_block(params->img, params->buf, params->w, params->h,
                             block);
        }
    }
    if (view->nthreads) thpool_wait(view->pool);
    uint64_t timer_stop = get_timeus();
    view->stats_ct += timer_stop - timer_start;
    view->stats_sps = 1000000 * (double)(view->cs * view->w * view->h +
                                         view->cb * view->bs * view->bs) /
                      (double)view->stats_ct;
    yg_update_texture(view->tid, view->img, view->w, view->h, view->nc);
    if (view->cb == view->nb) {
        view->cb = 0;
        view->cs++;
    }
}

// uiloop
void
ui_loop(const char* filename, const char* imfilename, yo_scene* scene,
        yt_scene* trace_scene, int w, int h, int ns, int bs, bool camera_lights,
        int stype, int rtype, float amb) {
    // view data
    view_params* view =
        init_view_params(filename, imfilename, scene, trace_scene, w, h, ns, bs,
                         camera_lights, stype, rtype, amb);

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

    thpool_destroy(view->pool);

    free(view);
}

#endif

bool
intersect_ray(void* ctx, const float ray_o[3], const float ray_d[3],
              float ray_min, float ray_max, int ray_mask, float* ray_t,
              int* sid, int* eid, float* euv) {
    yb_bvh* bvh = (yb_bvh*)ctx;
    return yb_intersect_bvh(bvh, ray_o, ray_d, ray_min, ray_max, ray_mask,
                            ray_t, sid, eid, euv);
}

bool
hit_ray(void* ctx, const float ray_o[3], const float ray_d[3], float ray_min,
        float ray_max, int ray_mask) {
    yb_bvh* bvh = (yb_bvh*)ctx;
    return yb_hit_bvh(bvh, ray_o, ray_d, ray_min, ray_max, ray_mask);
}

yb_bvh*
make_bvh(const yo_scene* scene) {
    // build bvh for rendering
    yb_bvh** shape_bvhs = (yb_bvh**)calloc(scene->nshapes, sizeof(yb_bvh*));
    int nshape_bvhs = 0;
    for (int sid = 0; sid < scene->nshapes; sid++) {
        yo_shape* shape = &scene->shapes[sid];
        shape_bvhs[nshape_bvhs] = yb_make_shape_bvh(
            shape->nelems, shape->elem, shape->etype, shape->nverts, shape->pos,
            shape->radius, 0, false);
        nshape_bvhs++;
    }
    // make xforms
    float** xforms = (float**)calloc(scene->nshapes, sizeof(float*));
    for (int i = 0; i < scene->nshapes; i++) {
        xforms[i] = scene->shapes[i].xform;
    }
    yb_bvh* scene_bvh =
        yb_make_scene_bvh(nshape_bvhs, shape_bvhs, xforms, 0, 0);
    free(xforms);
    return scene_bvh;
}

yt_scene*
init_trace_cb(const yo_scene* scene, yb_bvh* bvh, int camera) {
    yt_scene* trace_scene =
        yt_make_scene(scene->nshapes, scene->nmaterials, scene->ntextures);

    yt_set_intersection(trace_scene, bvh, intersect_ray, hit_ray);

    yo_camera* cam = &scene->cameras[camera];
    ym_mat4f camera_xform = ym_lookat_xform4f(
        *(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to, *(ym_vec3f*)cam->up);
    yt_set_camera(trace_scene, camera_xform.m, cam->width, cam->height,
                  cam->aperture,
                  ym_dist3f(*(ym_vec3f*)cam->from, *(ym_vec3f*)cam->to));

    if (scene->nenvs) {
        yo_env* env = &scene->envs[0];
        ym_mat4f env_xform = ym_lookat_xform4f(
            *(ym_vec3f*)env->from, *(ym_vec3f*)env->to, *(ym_vec3f*)env->up);
        yo_material* mat = scene->materials + env->matid;
        yt_set_env(trace_scene, env_xform.m, mat->ke, mat->ke_txtid);
    }

    for (int sid = 0; sid < scene->nshapes; sid++) {
        yo_shape* shape = &scene->shapes[sid];
        yt_set_shape(trace_scene, sid, shape->matid, shape->xform,
                     shape->nelems, shape->elem, shape->etype, shape->nverts,
                     shape->pos, shape->norm, shape->texcoord, 0, 0);
    }

    for (int mid = 0; mid < scene->nmaterials; mid++) {
        yo_material* mat = &scene->materials[mid];
        yt_set_material(trace_scene, mid, mat->ke, mat->kd, mat->ks,
                        yt_specular_exponent_to_roughness(mat->ns), 0,
                        mat->ke_txtid, mat->kd_txtid, mat->ks_txtid,
                        mat->ns_txtid, false);
    }

    for (int tid = 0; tid < scene->ntextures; tid++) {
        yo_texture* txt = &scene->textures[tid];
        yt_set_texture(trace_scene, tid, txt->pixels, txt->width, txt->height,
                       txt->ncomp);
    }

    yt_init_lights(trace_scene);

    return trace_scene;
}

void
render(yt_scene* trace_scene, const char* filename, const char* imfilename,
       int w, int h, int ns, int bs) {
    float* pixels = (float*)calloc(w * h * 4, sizeof(float));
    float* buf = (float*)calloc(w * h * 4, sizeof(float));
    int nblocks;
    int* blocks = make_image_blocks(w, h, bs, &nblocks);
    render_block_params* params =
        (render_block_params*)calloc(nblocks, sizeof(render_block_params));
    threadpool pool = thpool_init(4);
    printf("tracing %s to %s\n", filename, imfilename);
    printf("rendering ...");
    fflush(stdout);
    for (int s = 0; s < ns; s++) {
        printf("\rrendering sample %d/%d", s + 1, ns);
        fflush(stdout);
        for (int i = 0; i < nblocks; i++) {
            int* b = blocks + i * 4;
            params[i] = (render_block_params){
                pixels,     buf, w, h, ns, { b[0], b[1], b[2], b[3], s, s + 1 },
                trace_scene
            };
            thpool_add_work(pool, render_block_async, &params[i]);
        }
        thpool_wait(pool);
    }
    printf("\rrendering done\n");
    fflush(stdout);
    thpool_destroy(pool);
    save_image(imfilename, pixels, w, h, 4);
    free(buf);
    free(params);
    free(blocks);
    free(pixels);
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

    // ensure radius is necessary
    for (int s = 0; s < scene->nshapes; s++) {
        yo_shape* shape = &scene->shapes[s];
        if (shape->etype != yo_etype_point && shape->etype != yo_etype_line)
            continue;
        if (shape->radius) continue;
        shape->radius = (float*)calloc(shape->nverts, sizeof(float));
        for (int i = 0; i < shape->nverts; i++) shape->radius[i] = 0.001f;
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
    const char* rtype_names[] = { "default", "uniform", "stratified", "cmjs",
                                  0 };
    const char* stype_names[] = { "default", "eye", "direct", "path", 0 };
    // params
    yc_parser* parser = yc_init_parser(argc, argv, "trace meshes");
    int rtype =
        yc_parse_opte(parser, "--random", 0, "random type", 0, rtype_names);
    int stype = yc_parse_opte(parser, "--integrator", "-i", "integrator type",
                              0, stype_names);
    float amb = yc_parse_optf(parser, "--ambient", 0, "ambient factor", 0);
    bool triangulate =
        yc_parse_optb(parser, "--triangulate", 0, "triangulate input", false);
    bool camera_lights = yc_parse_optb(parser, "--camera_lights", "-c",
                                       "enable camera lights", false);
    int camera = yc_parse_opti(parser, "--camera", "-C", "camera", 0);
#ifndef YA_NOGL
    bool no_ui = yc_parse_optb(parser, "--no-ui", 0, "runs offline", false);
#endif
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

    // loading scene
    yo_scene* scene = load_scene(filename, triangulate);
    scene->cameras[camera].width = aspect * scene->cameras[camera].height;

    // preparing raytracer
    yb_bvh* scene_bvh = make_bvh(scene);
    yt_scene* trace_scene = init_trace_cb(scene, scene_bvh, camera);
    yt_set_rendering_params(trace_scene,
                            (camera_lights) ? yt_stype_eyelight : stype, rtype,
                            (float[3]){ amb, amb, amb });

#ifndef YA_NOGL
    // lunching renderer
    if (no_ui) {
        render(trace_scene, filename, imfilename, round(res * aspect), res,
               samples, block_size);
    } else {
        ui_loop(filename, imfilename, scene, trace_scene, round(res * aspect),
                res, samples, block_size, camera_lights, stype, rtype, amb);
    }
#else
    render(trace_scene, filename, imfilename, round(res * aspect), res, samples,
           block_size);
#endif

    // done
    // TODO: free bvh
    yt_free_scene(trace_scene);
    yo_free_scene(scene);
    return EXIT_SUCCESS;
}
