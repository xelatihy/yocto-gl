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

#include "../yocto/yocto_gl.h"
using namespace ygl;

// Generic image that contains either an HDR or an LDR image, giving access
// to both. This is helpful when writing viewers or generic image
// manipulation code
struct gimage {
    /// image path
    string filename;
    /// HDR image content
    image4f hdr;
    /// LDR image content
    image4b ldr;

    /// image width
    int width() const {
        if (!hdr.empty()) return hdr.width();
        if (!ldr.empty()) return ldr.width();
        return 0;
    }

    /// image height
    int height() const {
        if (!hdr.empty()) return hdr.height();
        if (!ldr.empty()) return ldr.height();
        return 0;
    }
};

// Loads a generic image
inline gimage load_gimage(const string& filename) {
    auto img = gimage();
    img.filename = filename;
    if (is_hdr_filename(filename)) {
        img.hdr = load_image4f(filename);
    } else {
        img.ldr = load_image4b(filename);
    }
    if (img.hdr.empty() && img.ldr.empty()) {
        throw runtime_error("cannot load image " + img.filename);
    }
    return img;
}

struct app_state {
    vector<gimage*> imgs;
    int cur_img = 0;

    gl_stdimage_program gl_prog = {};
    unordered_map<gimage*, gl_texture> gl_txt = {};
    gl_stdimage_params params;

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);
    auto img = app->imgs[app->cur_img];
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    gl_set_viewport(framebuffer_size);
    app->params.win_size = window_size;
    draw_image(app->gl_prog, app->gl_txt.at(img), app->params);

    if (begin_widgets(win, "yimview")) {
        draw_label_widget(win, "filename", img->filename);
        draw_label_widget(win, "size", "{} x {}", img->width(), img->height());
        draw_imageview_widgets(win, "", app->params, !img->hdr.empty());
        draw_imageinspect_widgets(
            win, "", img->hdr, img->ldr, get_mouse_posf(win), app->params);
    }
    end_widgets(win);

    swap_buffers(win);
}

void run_ui(app_state* app) {
    // window
    auto win = make_window(
        app->imgs[0]->width(), app->imgs[0]->height(), "yimview", app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // window values
    int mouse_button = 0;
    vec2f mouse_pos, mouse_last;

    // init widgets
    init_widgets(win);

    // load textures
    app->gl_prog = make_stdimage_program();
    for (auto img : app->imgs) {
        if (!img->hdr.empty()) {
            app->gl_txt[img] = make_texture(img->hdr, false, false, true);
        } else if (!img->ldr.empty()) {
            app->gl_txt[img] = make_texture(img->ldr, false, false, true);
        }
    }

    while (!should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = get_mouse_posf(win);
        mouse_button = get_mouse_button(win);

        auto img = app->imgs[app->cur_img];
        set_window_title(win, format("yimview | {} | {}x{}", img->filename,
                                  img->width(), img->height()));

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !get_widget_active(win)) {
            switch (mouse_button) {
                case 1: app->params.offset += mouse_pos - mouse_last; break;
                case 2:
                    app->params.zoom *=
                        powf(2, (mouse_pos[0] - mouse_last[0]) * 0.001f);
                    break;
                default: break;
            }
        }

        // draw
        draw(win);

        // event hadling
        wait_events(win);
    }

    clear_window(win);
    delete win;
}

int main(int argc, char* argv[]) {
    auto app = new app_state();

    // command line params
    auto parser = make_parser(argc, argv, "yimview", "view images");
    app->params.exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->params.gamma =
        parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->params.filmic =
        parse_flag(parser, "--filmic", "-F", "hdr image filmic");
    auto filenames =
        parse_args(parser, "image", "image filename", vector<string>{});
    // check parsing
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // loading images
    for (auto filename : filenames) {
        log_info("loading {}", filename);
        app->imgs.push_back(new gimage(load_gimage(filename)));
    }

    // run ui
    run_ui(app);

    // done
    delete app;
    return 0;
}
