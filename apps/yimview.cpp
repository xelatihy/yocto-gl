//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "../yocto/yocto_glutils.h"
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_utils.h"

#include <unordered_map>

// Generic image that contains either an HDR or an LDR image, giving access
// to both. This is helpful when writing viewers or generic image
// manipulation code
struct gimage {
    std::string filename;  // path

    int width = 0;                // width
    int height = 0;               // height
    std::vector<ygl::vec4f> hdr;  // hdr
    std::vector<ygl::vec4b> ldr;  // ldr

    // image adjustment
    bool updated = true;
    ygl::tonemap_type tonemap = ygl::tonemap_type::gamma;
    float exposure = 0;
    std::vector<ygl::vec4b> img;

    // opengl texture
    bool gl_updated = true;
    ygl::gltexture gl_txt;
};

// Loads a generic image
inline gimage* load_gimage(const std::string& filename) {
    auto img = new gimage();
    img->filename = filename;
    if (ygl::is_hdr_filename(filename)) {
        img->hdr = ygl::load_image4f(filename, img->width, img->height);
    } else {
        img->ldr = ygl::load_image4b(filename, img->width, img->height);
    }
    if (img->hdr.empty() && img->ldr.empty()) {
        throw std::runtime_error("cannot load image " + img->filename);
    }
    return img;
}

struct app_state {
    std::vector<gimage*> imgs;
    int cur_img = 0;

    ygl::glimage_program gl_prog = {};

    ygl::vec2f offset = {0, 0};
    float zoom = 1;
    ygl::vec4f background = {0, 0, 0, 0};

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

void draw_widgets(ygl::glwindow* win, app_state* app) {
    static auto tonemap_names = std::map<ygl::tonemap_type, std::string>{
        {ygl::tonemap_type::linear, "linear"},
        {ygl::tonemap_type::gamma, "gamma"},
        {ygl::tonemap_type::srgb, "srgb"},
        {ygl::tonemap_type::filmic1, "filmic1"},
        {ygl::tonemap_type::filmic2, "filmic2"},
        {ygl::tonemap_type::filmic3, "filmic3"},
    };

    if (ygl::begin_imgui_frame(win, "yimview")) {
        auto img = app->imgs[app->cur_img];
        ygl::draw_imgui_label(win, "filename", img->filename);
        ygl::draw_imgui_label(
            win, "size", ygl::format("{} x {}", img->width, img->height));
        ygl::draw_imgui_dragbox(win, "offset", app->offset, -4096, 4096);
        ygl::draw_imgui_dragbox(win, "zoom", app->zoom, 0.01, 10);
        ygl::draw_imgui_colorbox(win, "background", app->background);
        auto edited = 0;
        edited += ygl::draw_imgui_combobox(
            win, "tonemap", img->tonemap, tonemap_names);
        edited +=
            ygl::draw_imgui_dragbox(win, "exposure", img->exposure, -5, 5);
        img->updated = edited;
        auto ij = ygl::get_glimage_coords(
            get_glmouse_posf(win), app->offset, app->zoom);
        ygl::draw_imgui_dragbox(win, "mouse", ij);
        if (ij.x >= 0 && ij.x < img->width && ij.y >= 0 && ij.y < img->height) {
            if (!img->ldr.empty())
                ygl::draw_imgui_colorbox(
                    win, "pixel b", img->ldr.at(ij.x + ij.y * img->width));
            if (!img->hdr.empty())
                ygl::draw_imgui_colorbox(
                    win, "pixel f", img->hdr.at(ij.x + ij.y * img->width));
        } else {
            auto zero4b_ = ygl::zero4b;
            auto zero4f_ = ygl::zero4f;
            if (!img->ldr.empty())
                ygl::draw_imgui_colorbox(win, "pixel b", zero4b_);
            if (!img->hdr.empty())
                ygl::draw_imgui_colorbox(win, "pixel f", zero4f_);
        }
    }
    ygl::end_imgui_frame(win);
}

void draw(ygl::glwindow* win, app_state* app) {
    auto img = app->imgs[app->cur_img];
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glframebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::clear_glbuffers(app->background);
    ygl::draw_glimage(
        app->gl_prog, img->gl_txt, window_size, app->offset, app->zoom);
    draw_widgets(win, app);
    ygl::swap_glwindow_buffers(win);
}

void refresh(ygl::glwindow* win) {
    return draw(win, (app_state*)ygl::get_glwindow_user_pointer(win));
}

void run_ui(app_state* app) {
    // window
    auto win = ygl::make_glwindow(
        app->imgs[0]->width, app->imgs[0]->height, "yimview", app);
    ygl::set_glwindow_callbacks(win, nullptr, nullptr, refresh);

    // window values
    int mouse_button = 0;
    ygl::vec2f mouse_pos, mouse_last;

    // init widgets
    ygl::init_imgui(win);

    // load textures
    app->gl_prog = ygl::make_glimage_program();

    while (!should_glwindow_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = ygl::get_glmouse_posf(win);
        mouse_button = ygl::get_glmouse_button(win);

        auto img = app->imgs[app->cur_img];
        ygl::set_glwindow_title(
            win, ygl::format("yimview | {} | {}x{}", img->filename, img->width,
                     img->height));

        // handle mouse
        auto alt_down = get_glalt_key(win);
        if (mouse_button && alt_down && mouse_pos != mouse_last &&
            !ygl::get_imgui_active(win)) {
            switch (mouse_button) {
                case 1: app->offset += mouse_pos - mouse_last; break;
                case 2:
                    app->zoom *= powf(2, (mouse_pos.x - mouse_last.x) * 0.001f);
                    break;
                default: break;
            }
        }

        // update texture
        if (img->updated) {
            if (!img->hdr.empty()) {
                img->img =
                    ygl::tonemap_image(img->hdr, img->tonemap, img->exposure);
            } else {
                img->img = img->ldr;
            }
            img->gl_updated = true;
            img->updated = false;
        }
        if (img->gl_updated) {
            update_gltexture(img->gl_txt, img->width, img->height, img->img,
                false, false, false);
            img->gl_updated = false;
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
