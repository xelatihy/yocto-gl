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
    std::string name;
    std::string filename;  // path

    int width = 0;                // width
    int height = 0;               // height
    std::vector<ygl::vec4f> hdr;  // hdr
    std::vector<ygl::vec4b> ldr;  // ldr

    // min/max values
    ygl::vec4f hdr_min = ygl::zero4f, hdr_max = ygl::zero4f;
    ygl::vec4b ldr_min = ygl::zero4b, ldr_max = ygl::zero4b;
    float hdr_max_lum = 0;

    // image adjustment
    bool updated = true;
    ygl::tonemap_type tonemap = ygl::tonemap_type::gamma;
    float exposure = 0;
    std::vector<ygl::vec4b> img;

    // opengl texture
    bool gl_updated = true;
    ygl::gltexture gl_txt;
};

struct app_state {
    std::vector<gimage*> imgs;
    gimage* img = nullptr;

    ygl::glimage_program gl_prog = {};
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = false;
    ygl::vec4b background = {222, 222, 222, 0};

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

// compute min/max
void update_minmax(gimage* img) {
    if (!img->hdr.empty()) {
        img->hdr_min = {ygl::flt_max, ygl::flt_max, ygl::flt_max, ygl::flt_max};
        img->hdr_max = {ygl::flt_min, ygl::flt_min, ygl::flt_min, ygl::flt_min};
        img->hdr_max_lum = ygl::flt_min;
        for (auto& p : img->hdr) {
            img->hdr_min = {ygl::min(img->hdr_min.x, p.x),
                ygl::min(img->hdr_min.y, p.y), ygl::min(img->hdr_min.z, p.z),
                ygl::min(img->hdr_min.w, p.w)};
            img->hdr_max = {ygl::max(img->hdr_max.x, p.x),
                ygl::max(img->hdr_max.y, p.y), ygl::max(img->hdr_max.z, p.z),
                ygl::max(img->hdr_max.w, p.w)};
            img->hdr_max_lum = ygl::max(img->hdr_max_lum, luminance(p));
        }
    }
    if (!img->ldr.empty()) {
        img->ldr_min = {255, 255, 255, 255};
        img->ldr_max = {0, 0, 0, 0};
        for (auto& p : img->ldr) {
            img->ldr_min = {(ygl::byte)ygl::min(img->ldr_min.x, p.x),
                (ygl::byte)ygl::min(img->ldr_min.y, p.y),
                (ygl::byte)ygl::min(img->ldr_min.z, p.z),
                (ygl::byte)ygl::min(img->ldr_min.w, p.w)};
            img->ldr_max = {(ygl::byte)ygl::max(img->ldr_max.x, p.x),
                (ygl::byte)ygl::max(img->ldr_max.y, p.y),
                (ygl::byte)ygl::max(img->ldr_max.z, p.z),
                (ygl::byte)ygl::max(img->ldr_max.w, p.w)};
        }
    }
}

// Loads a generic image
gimage* load_gimage(const std::string& filename) {
    auto img = new gimage();
    img->filename = filename;
    img->name = ygl::path_filename(filename);
    if (ygl::is_hdr_filename(filename)) {
        img->hdr = ygl::load_image4f(filename, img->width, img->height);
    } else {
        img->ldr = ygl::load_image4b(filename, img->width, img->height);
    }
    if (img->hdr.empty() && img->ldr.empty()) {
        throw std::runtime_error("cannot load image " + img->filename);
    }
    update_minmax(img);
    return img;
}

gimage* diff_gimage(gimage* a, gimage* b, bool color) {
    if (a->width != b->width || a->height != b->height) return nullptr;
    if (a->ldr.empty() && !b->ldr.empty()) return nullptr;
    if (a->hdr.empty() && !b->hdr.empty()) return nullptr;
    auto d = new gimage();
    d->name = "diff " + a->name + " " + b->name;
    d->filename = "";
    d->width = a->width;
    d->height = a->height;
    if (!a->hdr.empty()) {
        d->hdr.resize(a->hdr.size());
        if (color) {
            for (auto i = 0; i < a->hdr.size(); i++) {
                d->hdr[i] = {std::abs(a->hdr[i].x - b->hdr[i].x),
                    std::abs(a->hdr[i].y - b->hdr[i].y),
                    std::abs(a->hdr[i].z - b->hdr[i].z),
                    std::max(a->hdr[i].w, b->hdr[i].w)};
            }
        } else {
            for (auto i = 0; i < a->hdr.size(); i++) {
                auto ld = std::abs(
                    ygl::luminance(a->hdr[i]) - ygl::luminance(b->hdr[i]));
                d->hdr[i] = {ld, ld, ld, std::max(a->hdr[i].w, b->hdr[i].w)};
            }
        }
    }
    if (!a->ldr.empty()) {
        d->ldr.resize(a->ldr.size());
        for (auto i = 0; i < a->ldr.size(); i++)
            d->ldr[i] = {(ygl::byte)std::abs(a->ldr[i].x - b->ldr[i].x),
                (ygl::byte)std::abs(a->ldr[i].y - b->ldr[i].y),
                (ygl::byte)std::abs(a->ldr[i].z - b->ldr[i].z),
                std::max(a->ldr[i].w, b->ldr[i].w)};
    }
    update_minmax(d);
    return d;
}

void update_display_image(gimage* img) {
    if (!img->hdr.empty()) {
        img->img = ygl::tonemap_image(img->hdr, img->tonemap, img->exposure);
    } else {
        img->img = img->ldr;
    }
    img->updated = false;
}

void draw_glwidgets(ygl::glwindow* win, app_state* app) {
    static auto tonemap_names = std::map<ygl::tonemap_type, std::string>{
        {ygl::tonemap_type::linear, "linear"},
        {ygl::tonemap_type::gamma, "gamma"},
        {ygl::tonemap_type::srgb, "srgb"},
        {ygl::tonemap_type::filmic1, "filmic1"},
        {ygl::tonemap_type::filmic2, "filmic2"},
        {ygl::tonemap_type::filmic3, "filmic3"},
    };

    if (ygl::begin_glwidgets_frame(win, "yimview")) {
        ygl::draw_glwidgets_combobox(win, "image", app->img, app->imgs);
        ygl::draw_glwidgets_label(win, "filename", app->img->filename);
        ygl::draw_glwidgets_label(win, "size",
            ygl::format("{} x {}", app->img->width, app->img->height));
        auto edited = 0;
        edited += ygl::draw_glwidgets_combobox(
            win, "tonemap", app->img->tonemap, tonemap_names);
        edited += ygl::draw_glwidgets_dragbox(
            win, "exposure", app->img->exposure, -5, 5);
        if (edited) app->img->updated = true;
        auto zoom = app->imframe.x.x;
        if (ygl::draw_glwidgets_dragbox(win, "zoom", zoom, 0.1, 10))
            app->imframe.x.x = app->imframe.y.y = zoom;
        ygl::draw_glwidgets_checkbox(win, "zoom to fit", app->zoom_to_fit);
        ygl::draw_glwidgets_colorbox(win, "background", app->background);
        ygl::draw_glwidgets_separator(win);
        auto ij = ygl::get_glimage_coords(get_glwidnow_mouse_posf(win),
            app->imframe, {app->img->width, app->img->height});
        ygl::draw_glwidgets_dragbox(win, "mouse", ij);
        if (ij.x >= 0 && ij.x < app->img->width && ij.y >= 0 &&
            ij.y < app->img->height) {
            if (!app->img->ldr.empty()) {
                ygl::draw_glwidgets_colorbox(win, "pixel",
                    app->img->ldr.at(ij.x + ij.y * app->img->width));
            }
            if (!app->img->hdr.empty()) {
                ygl::draw_glwidgets_colorbox(win, "pixel",
                    app->img->hdr.at(ij.x + ij.y * app->img->width));
            }
        } else {
            auto zero4b_ = ygl::zero4b;
            auto zero4f_ = ygl::zero4f;
            if (!app->img->ldr.empty())
                ygl::draw_glwidgets_colorbox(win, "pixel", zero4b_);
            if (!app->img->hdr.empty())
                ygl::draw_glwidgets_colorbox(win, "pixel", zero4f_);
        }
        if (!app->img->ldr.empty()) {
            ygl::draw_glwidgets_colorbox(win, "min", app->img->ldr_min);
            ygl::draw_glwidgets_colorbox(win, "max", app->img->ldr_max);
        }
        if (!app->img->hdr.empty()) {
            ygl::draw_glwidgets_colorbox(win, "min", app->img->hdr_min);
            ygl::draw_glwidgets_colorbox(win, "max", app->img->hdr_max);
            float max_exposure = log(app->img->hdr_max_lum + 0.00001f) / log(2);
            ygl::draw_glwidgets_dragbox(win, "max l", app->img->hdr_max_lum);
            ygl::draw_glwidgets_dragbox(win, "max e", max_exposure);
        }
    }
    ygl::end_glwidgets_frame(win);
}

void draw(ygl::glwindow* win, app_state* app) {
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glwindow_framebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::clear_glbuffers(app->background);
    ygl::draw_glimage(
        app->gl_prog, app->img->gl_txt, window_size, app->imframe);
    draw_glwidgets(win, app);
    ygl::swap_glwindow_buffers(win);
}

void refresh(ygl::glwindow* win) {
    auto app = (app_state*)ygl::get_glwindow_user_pointer(win);
    ygl::center_glimage(app->imframe, {app->img->width, app->img->height},
        ygl::get_glwindow_size(win), app->zoom_to_fit);
    draw(win, app);
}

void run_ui(app_state* app) {
    // window
    auto win_width = app->imgs[0]->width + ygl::default_glwidgets_width;
    auto win_height = ygl::clamp(app->imgs[0]->height, 512, 1024);
    auto win = ygl::make_glwindow(win_width, win_height, "yimview", app);
    ygl::set_glwindow_callbacks(win, nullptr, nullptr, refresh);
    ygl::center_glimage(app->imframe,
        {app->imgs[0]->width, app->imgs[0]->height},
        ygl::get_glwindow_size(win), app->zoom_to_fit);

    // window values
    int mouse_button = 0;
    ygl::vec2f mouse_pos, mouse_last;

    // init widgets
    ygl::init_glwidgets(win);

    // load textures
    app->gl_prog = ygl::make_glimage_program();

    while (!should_glwindow_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = ygl::get_glwidnow_mouse_posf(win);
        mouse_button = ygl::get_glwindow_mouse_button(win);

        ygl::set_glwindow_title(
            win, ygl::format("yimview | {} | {}x{}", app->img->filename,
                     app->img->width, app->img->height));

        // handle mouse
        auto alt_down = get_glwindow_alt_key(win);
        if (mouse_button && alt_down && mouse_pos != mouse_last &&
            !ygl::get_glwidgets_active(win)) {
            switch (mouse_button) {
                case 1: app->imframe.o += mouse_pos - mouse_last; break;
                case 2: {
                    auto zoom = powf(2, (mouse_pos.x - mouse_last.x) * 0.001f);
                    app->imframe = app->imframe *
                                   ygl::frame2f{{zoom, 0}, {0, zoom}, {0, 0}};
                } break;
                default: break;
            }
        }
        // auto center
        ygl::center_glimage(app->imframe, {app->img->width, app->img->height},
            ygl::get_glwindow_size(win), app->zoom_to_fit);

        // update texture
        if (app->img->updated) {
            update_display_image(app->img);
            app->img->gl_updated = true;
        }
        if (app->img->gl_updated) {
            update_gltexture(app->img->gl_txt, app->img->width,
                app->img->height, app->img->img, false, false);
            app->img->gl_updated = false;
        }

        // draw
        draw(win, app);

        // event hadling
        if (ygl::get_glwindow_mouse_button(win) ||
            ygl::get_glwidgets_active(win)) {
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
    auto diff = ygl::parse_flag(parser, "--diff", "-d", "compute diff images");
    auto lum_diff = ygl::parse_flag(
        parser, "--luminance-diff", "-D", "compute luminance diffs");
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
    app->img = app->imgs.at(0);
    if (diff) {
        auto num = app->imgs.size();
        for (auto i = 0; i < num; i++) {
            for (auto j = i + 1; j < num; j++) {
                ygl::log_info("diffing {} {}", app->imgs[i]->filename,
                    app->imgs[j]->filename);
                auto diff = diff_gimage(app->imgs[i], app->imgs[j], !lum_diff);
                if (diff) app->imgs.push_back(diff);
            }
        }
    }

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
