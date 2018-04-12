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
        if (!hdr.pixels.empty()) return hdr.width;
        if (!ldr.pixels.empty()) return ldr.width;
        return 0;
    }

    /// image height
    int height() const {
        if (!hdr.pixels.empty()) return hdr.height;
        if (!ldr.pixels.empty()) return ldr.height;
        return 0;
    }
};

// Loads a generic image
inline gimage* load_gimage(const std::string& filename) {
    auto img = new gimage();
    img->filename = filename;
    if (ygl::is_hdr_filename(filename)) {
        img->hdr = ygl::load_image4f(filename);
    } else {
        img->ldr = ygl::load_image4b(filename);
    }
    if (img->hdr.pixels.empty() && img->ldr.pixels.empty()) {
        throw std::runtime_error("cannot load image " + img->filename);
    }
    return img;
}

struct app_state {
    std::vector<gimage*> imgs;
    int cur_img = 0;

    ygl::glimage_program gl_prog = {};
    std::unordered_map<gimage*, ygl::gltexture> gl_txt = {};
    ygl::glimage_params params;
    float exposure = 1.0f;

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

void draw(ygl::glwindow* win, app_state* app) {
    auto img = app->imgs[app->cur_img];
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glframebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::draw_glimage(
        app->gl_prog, app->gl_txt.at(img), window_size, app->params);

    if (ygl::begin_imgui_frame(win, "yimview")) {
        ygl::draw_imgui_label(win, "filename", img->filename);
        ygl::draw_imgui_label(
            win, "size", "{} x {}", img->width(), img->height());
        ygl::draw_imgui_stdimage_inspector(win, "", app->params);
        ygl::draw_imgui_image_inspector(
            win, "", img->hdr, img->ldr, get_glmouse_posf(win), app->params);
    }
    ygl::end_imgui_frame(win);

    ygl::swap_glwindow_buffers(win);
}

void run_ui(app_state* app) {
    // window
    auto win = ygl::make_glwindow(
        app->imgs[0]->width(), app->imgs[0]->height(), "yimview");
    ygl::set_glwindow_callbacks(
        win, nullptr, nullptr, [app, win]() { draw(win, app); });

    // window values
    int mouse_button = 0;
    ygl::vec2f mouse_pos, mouse_last;

    // init widgets
    ygl::init_imgui(win);

    // load textures
    app->gl_prog = ygl::make_glimage_program();
    for (auto img : app->imgs) {
        if (!img->hdr.pixels.empty()) {
            update_gltexture(app->gl_txt[img], img->hdr, false, false, true);
        } else if (!img->ldr.pixels.empty()) {
            update_gltexture(app->gl_txt[img], img->ldr, false, false, true);
        }
    }

    while (!should_glwindow_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = ygl::get_glmouse_posf(win);
        mouse_button = ygl::get_glmouse_button(win);

        auto img = app->imgs[app->cur_img];
        ygl::set_glwindow_title(
            win, ygl::format("yimview | {} | {}x{}", img->filename,
                     img->width(), img->height()));

        // handle mouse
        auto alt_down = get_glalt_key(win);
        if (mouse_button && alt_down && mouse_pos != mouse_last &&
            !ygl::get_imgui_active(win)) {
            switch (mouse_button) {
                case 1: app->params.offset += mouse_pos - mouse_last; break;
                case 2:
                    app->params.zoom *=
                        powf(2, (mouse_pos.x - mouse_last.x) * 0.001f);
                    break;
                default: break;
            }
        }

        // draw
        draw(win, app);

        // event hadling
        if (ygl::get_glmouse_button(win) || ygl::get_imgui_active(win)) {
            ygl::poll_glwindow_events(win);
        } else {
            ygl::wait_glwindow_events(win);
        }
    }

    // cleanup
    delete win;
}

int main(int argc, char* argv[]) {
    auto app = new app_state();

    // command line params
    auto parser = ygl::make_parser(argc, argv, "yimview", "view images");
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
        app->imgs.push_back(load_gimage(filename));
    }

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
