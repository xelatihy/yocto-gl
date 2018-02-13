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

// Generic image that contains either an HDR or an LDR image, giving access
// to both. This is helpful when writing viewers or generic image
// manipulation code
struct gimage {
    /// image path
    std::string filename;
    /// HDR image content
    ygl::image4f hdr;
    /// LDR image content
    ygl::image4b ldr;

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
inline gimage load_gimage(const std::string& filename) {
    auto img = gimage();
    img.filename = filename;
    if (ygl::is_hdr_filename(filename)) {
        img.hdr = ygl::load_image4f(filename);
    } else {
        img.ldr = ygl::load_image4b(filename);
    }
    if (img.hdr.empty() && img.ldr.empty()) {
        throw std::runtime_error("cannot load image " + img.filename);
    }
    return img;
}

struct app_state {
    std::vector<gimage*> imgs;
    int cur_img = 0;

    ygl::gl_stdimage_program gl_prog = {};
    std::unordered_map<gimage*, ygl::gl_texture> gl_txt = {};
    ygl::gl_stdimage_params params;

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

void draw(ygl::gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);
    auto img = app->imgs[app->cur_img];
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    ygl::gl_set_viewport(framebuffer_size);
    ygl::draw_image(
        app->gl_prog, app->gl_txt.at(img), window_size, app->params);

    if (ygl::begin_widgets(win, "yimview")) {
        ygl::draw_label_widget(win, "filename", img->filename);
        ygl::draw_label_widget(
            win, "size", "{} x {}", img->width(), img->height());
        ygl::draw_params_widgets(win, "", app->params);
        ygl::draw_imageinspect_widgets(
            win, "", img->hdr, img->ldr, get_mouse_posf(win), app->params);
    }
    ygl::end_widgets(win);

    ygl::swap_buffers(win);
}

void run_ui(app_state* app) {
    // window
    auto win = ygl::make_window(
        app->imgs[0]->width(), app->imgs[0]->height(), "yimview", app);
    ygl::set_window_callbacks(win, nullptr, nullptr, draw);

    // window values
    int mouse_button = 0;
    ygl::vec2f mouse_pos, mouse_last;

    // init widgets
    ygl::init_widgets(win);

    // load textures
    app->gl_prog = ygl::make_stdimage_program();
    for (auto img : app->imgs) {
        if (!img->hdr.empty()) {
            app->gl_txt[img] = make_texture(img->hdr, false, false, true);
        } else if (!img->ldr.empty()) {
            app->gl_txt[img] = make_texture(img->ldr, false, false, true);
        }
    }

    while (!should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = ygl::get_mouse_posf(win);
        mouse_button = ygl::get_mouse_button(win);

        auto img = app->imgs[app->cur_img];
        ygl::set_window_title(
            win, ygl::format("yimview | {} | {}x{}", img->filename,
                     img->width(), img->height()));

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !ygl::get_widget_active(win)) {
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
        ygl::wait_events(win);
    }

    ygl::clear_window(win);
    delete win;
}

int main(int argc, char* argv[]) {
    auto app = new app_state();

    // command line params
    auto parser = ygl::make_parser(argc, argv, "yimview", "view images");
    app->params = ygl::parse_params(parser, "", app->params);
    auto filenames = ygl::parse_args(
        parser, "image", "image filename", std::vector<std::string>{});
    // check parsing
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // loading images
    for (auto filename : filenames) {
        ygl::log_info("loading {}", filename);
        app->imgs.push_back(new gimage(load_gimage(filename)));
    }

    // run ui
    run_ui(app);

    // done
    delete app;
    return 0;
}
