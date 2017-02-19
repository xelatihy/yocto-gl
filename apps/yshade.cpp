//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#include "yapp.h"

#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"

struct state {
    // params
    std::shared_ptr<yapp::params> pars = nullptr;

    // scene
    std::shared_ptr<yapp::scene> scene = nullptr;

    // shade state
    std::shared_ptr<yapp::shade_state> shst = nullptr;

    // widgets
    void* widget_ctx = nullptr;
};

const int hud_width = 256;

void save_screenshot(yglu::ui::window* win, const std::string& imfilename) {
    if (ycmd::get_extension(imfilename) != ".png") {
        printf("supports only png screenshots");
        return;
    }

    auto wh = yapp::int2{0, 0};
    auto pixels = yglu::ui::get_screenshot(win, wh);
    yimg::save_image(imfilename, wh[0], wh[1], 4, nullptr,
                     (unsigned char*)pixels.data());
}

void draw_scene(yglu::ui::window* win) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    auto window_size = yglu::ui::get_window_size(win);
    st->scene->cameras[pars->render_params.camera_id]->aspect =
        (float)window_size[0] / (float)window_size[1];
    yapp::shade_scene(st->scene, pars, st->shst);
}

void text_callback(yglu::ui::window* win, unsigned int key) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    switch (key) {
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '\\': pars->srgb = !pars->srgb; break;
        case 'w': pars->wireframe = !pars->wireframe; break;
        case 'e': pars->edges = !pars->edges; break;
        case 's': save_screenshot(win, pars->imfilename); break;
        case 'C':
            pars->render_params.camera_id =
                (pars->render_params.camera_id + 1) % st->scene->cameras.size();
            break;
        case 't': {
            for (auto shape : st->scene->shapes) {
                yshape::tesselate_stdshape(
                    (std::vector<yshape::int2>&)shape->lines,
                    (std::vector<yshape::int3>&)shape->triangles,
                    (std::vector<yshape::float3>&)shape->pos,
                    (std::vector<yshape::float3>&)shape->norm,
                    (std::vector<yshape::float2>&)shape->texcoord,
                    (std::vector<yshape::float3>&)shape->color, shape->radius);
            }
        } break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_widgets(yglu::ui::window* win) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    if (yglu::ui::begin_widgets(win)) {
        yglu::ui::dynamic_widget_layout(win, 1);
        yglu::ui::label_widget(win, pars->filenames[0]);
        yglu::ui::dynamic_widget_layout(win, 2);
        yglu::ui::int_widget(win, "camera", &pars->render_params.camera_id, 0,
                             (int)st->scene->cameras.size() - 1);
        yglu::ui::bool_widget(win, "wireframe", &pars->wireframe);
        yglu::ui::bool_widget(win, "edges", &pars->edges);
        yglu::ui::dynamic_widget_layout(win, 1);
        yglu::ui::float_widget(win, "exposure", &pars->exposure, -20, 20, 1);
        yglu::ui::float_widget(win, "gamma", &pars->gamma, 0.1, 5, 0.1);
        yglu::ui::bool_widget(win, "srgb output", &pars->srgb);
        if (yglu::ui::button_widget(win, "tesselate")) {
            for (auto shape : st->scene->shapes) {
                yshape::tesselate_stdshape(
                    (std::vector<yshape::int2>&)shape->lines,
                    (std::vector<yshape::int3>&)shape->triangles,
                    (std::vector<yshape::float3>&)shape->pos,
                    (std::vector<yshape::float3>&)shape->norm,
                    (std::vector<yshape::float2>&)shape->texcoord,
                    (std::vector<yshape::float3>&)shape->color, shape->radius);
            }
        }
    }
    yglu::ui::end_widgets(win);
}

void window_refresh_callback(yglu::ui::window* win) {
    draw_scene(win);
    draw_widgets(win);
    yglu::ui::swap_buffers(win);
}

void run_ui(const std::shared_ptr<state>& st) {
    auto pars = st->pars;

    // window
    auto win = yglu::ui::init_window(pars->width, pars->height, "yshade",
                                     pars->legacy_gl, st.get());
    yglu::ui::set_callbacks(win, text_callback, window_refresh_callback);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;

    // load textures
    st->shst = yapp::init_shade_state(st->scene, pars);

    // init widget
    yglu::ui::init_widgets(win);

    // loop
    while (!yglu::ui::should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = yglu::ui::get_mouse_posf(win);
        mouse_button = yglu::ui::get_mouse_button(win);

        yglu::ui::set_window_title(win, ("yshade | " + pars->filenames[0]));

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !yglu::ui::get_widget_active(win)) {
            auto dolly = 0.0f;
            auto pan = ym::zero2f;
            auto rotate = ym::zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }

            auto cam = st->scene->cameras[pars->render_params.camera_id];
            ym::turntable((ym::frame3f&)cam->frame, cam->focus, rotate, dolly,
                          pan);
        }

        // draw
        draw_scene(win);

        // make ui
        draw_widgets(win);

        // swap buffers
        yglu::ui::swap_buffers(win);

        // check for screenshot
        if (pars->no_ui) {
            save_screenshot(win, pars->imfilename);
            break;
        }

        // event hadling
        yglu::ui::wait_events(win);
    }

    yglu::ui::clear_widgets(win);
    yglu::ui::clear_window(win);
}

int main(int argc, char* argv[]) {
    // params
    auto pars = yapp::init_params("interactively view scenes", argc, argv,
                                  false, false, true, true);

    // init state
    auto st = std::make_shared<state>();
    st->pars = pars;

    // setting up rendering
    st->scene = yapp::load_scenes(pars->filenames, pars->scene_scale);

    // run ui
    run_ui(st);

    // done
    return 0;
}
