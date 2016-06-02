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

#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"

#ifndef STBI_INCLUDE_STB_IMAGE_H

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_STATIC

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#endif

#include "ext/stb_image.h"

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

#endif

struct view_img {
    ym_string filename;
    ym_vector<float> pixels;
    int w, h, nc;
    int tex_glid;
};
typedef struct view_img view_img;

struct view_params {
    int w, h, fw, fh;
    float zoom, ox, oy;
    float exposure, gamma;
    float* background;
    float backgrounds[16];

    bool hud_print;

    int mouse_button;
    float mouse_x, mouse_y;

    int nimgs;
    view_img* imgs;
    view_img* img;
};
typedef struct view_params view_params;

view_params* init_view_params(int nimgs, view_img* imgs) {
    view_params* view = new view_params();

    view->w = imgs[0].w;
    view->h = imgs[0].h;
    view->imgs = imgs;
    view->nimgs = nimgs;
    view->img = view->imgs;

    view->zoom = 1;
    view->gamma = 2.2f;
    view->exposure = 0;
    view->background = view->backgrounds;

    float backgrounds[16] = {0,    0,    0,    0, 0.18f, 0.18f, 0.18f, 0,
                             0.5f, 0.5f, 0.5f, 0, 1,     1,     1,     0};
    memcpy(view->backgrounds, backgrounds, sizeof(view->backgrounds));

    view->hud_print = false;

    return view;
}

// inspect string
#define strcatf(str, fmt, ...)                                                 \
    {                                                                          \
        char __buf[512];                                                       \
        sprintf(__buf, fmt, __VA_ARGS__);                                      \
        strcat(str, __buf);                                                    \
    }

void inspect_str(char* buf, view_params* view) {
    view_img* img = view->img;

    int ix = (int)round(((int)view->mouse_x - view->ox) / view->zoom);
    int iy = (int)round(((int)view->mouse_y - view->oy) / view->zoom);

    strcpy(buf, "");
    strcatf(buf, "img: %d x %d @ %d\n", img->w, img->h, img->nc);
    strcatf(buf, "mouse: %4d %4d\n", ix, iy);
    if (ix >= 0 && ix < img->w && iy >= 0 && iy < img->h) {
        float c[4] = {0, 0, 0, 1};
        int cb[4] = {0, 0, 0, 255};
        memcpy(c, img->pixels.data() + (iy * img->w + ix) * img->nc,
               sizeof(float) * img->nc);
        cb[0] = ym_clamp(
            (int)(255 * pow(2, view->exposure) * pow(c[0], 1 / view->gamma)), 0,
            255);
        cb[1] = ym_clamp(
            (int)(255 * pow(2, view->exposure) * pow(c[1], 1 / view->gamma)), 0,
            255);
        cb[2] = ym_clamp(
            (int)(255 * pow(2, view->exposure) * pow(c[2], 1 / view->gamma)), 0,
            255);
        cb[3] = ym_clamp((int)(255 * c[3]), 0, 255);
        strcatf(buf, "r: %10.4f %3d\n", c[0], cb[0]);
        strcatf(buf, "g: %10.4f %3d\n", c[1], cb[1]);
        strcatf(buf, "b: %10.4f %3d\n", c[2], cb[2]);
        strcatf(buf, "a: %10.4f %3d\n", c[3], cb[3]);
    } else {
        strcatf(buf, "r: %s\n", "-");
        strcatf(buf, "g: %s\n", "-");
        strcatf(buf, "b: %s\n", "-");
        strcatf(buf, "a: %s\n", "-");
    }
}

// text callback
void text_callback(GLFWwindow* window, unsigned int key) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    switch (key) {
        case ' ':
        case '.':
            view->img++;
            if (view->img >= view->imgs + view->nimgs) view->img = view->imgs;
            break;
        case ',':
            view->img--;
            if (view->img < view->imgs)
                view->img = view->imgs + view->nimgs - 1;
            break;
        case '-':
        case '_': view->zoom /= 2; break;
        case '+':
        case '=': view->zoom *= 2; break;
        case '[': view->exposure -= 1; break;
        case ']': view->exposure += 1; break;
        case '{': view->gamma -= 0.1f; break;
        case '}': view->gamma += 0.1f; break;
        case '1':
            view->exposure = 0;
            view->gamma = 1;
            break;
        case '2':
            view->exposure = 0;
            view->gamma = 2.2f;
            break;
        case 'z': view->zoom = 1; break;
        case 'b':
            view->background += 4;
            if (view->background - view->backgrounds >= 16)
                view->background = view->backgrounds;
            break;
        case 'h':
            // TODO: hud
            break;
        case 'p': view->hud_print = !view->hud_print; break;
        default: printf("unsupported key\n"); break;
    }
}

void window_size_callback(GLFWwindow* window, int w, int h) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);
    view->w = w;
    view->h = h;
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

    switch (view->mouse_button) {
        case 1:
            view->ox += (float)x - view->mouse_x;
            view->oy += (float)y - view->mouse_y;
            break;
        case 2:
            view->zoom *= powf(2, (view->mouse_x - (float)x) * 0.001f);
            break;
        default: break;
    }

    view->mouse_x = (float)x;
    view->mouse_y = (float)y;

    // draw inspector
    if (view->hud_print) {
        char buf[4096];
        inspect_str(buf, view);
        printf("\033[2J\033[1;1H");
        printf("%s", buf);
    }
}

void window_refresh_callback(GLFWwindow* window) {
    view_params* view = (view_params*)glfwGetWindowUserPointer(window);

    char title[4096];
    sprintf(title, "yview | %s | %dx%d", view->img->filename.c_str(),
            view->img->w, view->img->h);
    glfwSetWindowTitle(window, title);

    // grab selected image
    view_img* img = view->img;

    // begin frame
    glClearColor(view->background[0], view->background[1], view->background[2],
                 view->background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    yg_shade_image(img->tex_glid, img->w, img->h, view->w, view->h, view->ox,
                   view->oy, view->zoom, view->exposure, view->gamma);

    glfwSwapBuffers(window);
}

// uiloop
void ui_loop(int nimgs, view_img* imgs) {
    // view params
    view_params* view = init_view_params(nimgs, imgs);

    // window
    if (!glfwInit()) exit(EXIT_FAILURE);
    GLFWwindow* window = glfwCreateWindow(view->w, view->h, "imview", 0, 0);
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

    // load textures
    for (int i = 0; i < nimgs; i++) {
        imgs[i].tex_glid = yg_make_texture(imgs[i].pixels.data(), imgs[i].w,
                                           imgs[i].h, imgs[i].nc, true, false);
    }

    // ui loop
    while (!glfwWindowShouldClose(window)) {
        window_refresh_callback(window);
        glfwWaitEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    free(view);
}

ym_vector<view_img> load_images(const char** img_filenames) {
    ym_vector<view_img> imgs;
    for (const char** filename = img_filenames; *filename; filename++) {
        imgs.push_back(view_img());
        view_img* img = &imgs[imgs.size() - 1];
        img->filename = *filename;
        float* pixels =
            stbi_loadf(img->filename.c_str(), &img->w, &img->h, &img->nc, 0);
        img->pixels =
            ym_vector<float>(pixels, pixels + img->w * img->h * img->nc);
        free(pixels);
        if (img->pixels.empty()) {
            printf("cannot load image %s\n", img->filename.c_str());
            exit(1);
        }
        img->tex_glid = 0;
    }
    return imgs;
}

int main(int argc, const char** argv) {
    const char* filenames[4096];

    // command line params
    yc_parser* parser = yc_init_parser(argc, argv, "view images");
    yc_parse_argsa(parser, "image", "image filename", 0, true, filenames, 4096);
    yc_done_parser(parser);

    // loading images
    ym_vector<view_img> imgs = load_images(filenames);

    // start ui loop
    ui_loop((int)imgs.size(), imgs.data());

    // done
    return EXIT_SUCCESS;
}
