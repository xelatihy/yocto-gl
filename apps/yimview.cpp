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

#define YGL_OPENGL 1
#include "../yocto/yocto_gl.h"
using namespace ygl;

/// Generic image that contains either an HDR or an LDR image, giving access
/// to both. This is helpful when writing viewers or generic image
/// manipulation code
struct gimage {
    /// image path
    string filename;
    /// HDR image content
    image4f hdr;
    /// LDR image content
    image4b ldr;

    /// Check if the image is valid
    operator bool() const { return hdr || ldr; }

    /// image width
    int width() const {
        if (hdr) return hdr.width();
        if (ldr) return ldr.width();
        return 0;
    }

    /// image height
    int height() const {
        if (hdr) return hdr.height();
        if (ldr) return ldr.height();
        return 0;
    }

    /// access to pixel values
    vec4f& at4f(const vec2i& ij) { return hdr.at(ij); }
    /// access to pixel values
    const vec4f& at4f(const vec2i& ij) const { return hdr.at(ij); }
    /// access to pixel values
    vec4b& at4b(const vec2i& ij) { return ldr.at(ij); }
    /// access to pixel values
    const vec4b& at4b(const vec2i& ij) const { return ldr.at(ij); }

    /// guarded access to pixel values
    vec4f lookup4f(const vec2i& ij, const vec4f& def = zero4f) const {
        if (ij.x < 0 || ij.x >= width() || ij.y < 0 || ij.y > height())
            return def;
        if (hdr) return hdr.at(ij);
        if (ldr) return srgb_to_linear(ldr.at(ij));
        return def;
    }
    /// guarded access to pixel values
    vec4b lookup4b(const vec2i& ij, const vec4b& def = zero4b) const {
        if (ij.x < 0 || ij.x >= width() || ij.y < 0 || ij.y > height())
            return def;
        if (ldr) return ldr.at(ij);
        if (hdr) return linear_to_srgb(hdr.at(ij));
        return def;
    }
};

/// Loads a generic image
inline gimage load_gimage(const string& filename) {
    auto img = gimage();
    img.filename = filename;
    if (is_hdr_filename(filename)) {
        img.hdr = load_image4f(filename);
    } else {
        img.ldr = load_image4b(filename);
    }
    if (!img) { throw runtime_error("cannot load image " + img.filename); }
    return img;
}

struct app_state {
    vector<gimage*> imgs;
    int cur_img = 0;

    float exposure = 0;
    float gamma = 1;
    bool filmic = false;

    float zoom = 1;
    vec2f offset = vec2f();

    float background = 0;

    gl_stdimage_program gl_prog = {};
    unordered_map<gimage*, gl_texture> gl_txt = {};

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);
    auto img = app->imgs[app->cur_img];
    gl_clear_buffers({app->background, app->background, app->background, 0});
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    gl_set_viewport(framebuffer_size);
    draw_image(app->gl_prog, app->gl_txt.at(img), window_size, app->offset,
        app->zoom, app->exposure, app->gamma, app->filmic);

    auto mouse_pos = (vec2f)get_mouse_posf(win);
    if (begin_widgets(win, "yimview")) {
        draw_label_widget(win, "filename", img->filename);
        draw_label_widget(win, "w", img->width());
        draw_label_widget(win, "h", img->height());
        if (img->hdr) {
            draw_value_widget(win, "filmic", app->filmic);
            draw_value_widget(win, "exposure", app->exposure, -20, 20, 1);
            draw_value_widget(win, "gamma", app->gamma, 0.1, 5, 0.1);
        }
        draw_value_widget(win, "offset", app->offset, -100, 100, 1);
        draw_value_widget(win, "zoom", app->zoom, 0.01, 10, 0.1);
        auto xy = (mouse_pos - app->offset) / app->zoom;
        auto ij = vec2i{(int)round(xy[0]), (int)round(xy[1])};
        auto ph = img->lookup4f(ij);
        draw_label_widget(win, "r", ph.x);
        draw_label_widget(win, "g", ph.y);
        draw_label_widget(win, "b", ph.z);
        draw_label_widget(win, "a", ph.w);
        auto pl = img->lookup4b(ij);
        draw_label_widget(win, "r", pl.x);
        draw_label_widget(win, "g", pl.y);
        draw_label_widget(win, "b", pl.z);
        draw_label_widget(win, "a", pl.w);
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
        if (img->hdr) {
            app->gl_txt[img] = make_texture(img->hdr, false, false, true);
        } else if (img->ldr) {
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
                case 1: app->offset += mouse_pos - mouse_last; break;
                case 2:
                    app->zoom *=
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
    app->exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->gamma = parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->filmic = parse_flag(parser, "--filmic", "-F", "hdr image filmic");
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
