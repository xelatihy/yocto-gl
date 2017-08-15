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

#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_gui.h"
#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_utils.h"

using yu::logging::log_fatal;

struct yimage {
    // image path
    std::string filename;

    // original image data size
    int width() const {
        if (hdr) return hdr.width();
        if (ldr) return ldr.width();
        return 0;
    }
    int height() const {
        if (hdr) return hdr.height();
        if (ldr) return ldr.height();
        return 0;
    }

    // pixel data
    ym::image<ym::vec4f> hdr;
    ym::image<ym::vec4b> ldr;

    // opengl texture
    yglu::uint tex_glid = 0;
};

struct params {
    std::vector<std::string> filenames;
    std::vector<yimage> imgs;

    float exposure = 0;
    float gamma = 1;
    ym::tonemap_type tonemap = ym::tonemap_type::gamma;

    int cur_img = 0;
    int cur_background = 0;
    float zoom = 1;
    ym::vec2f offset = ym::vec2f();

    float background = 0;

    void* widget_ctx = nullptr;
};

bool load_images(params* pars) {
    for (auto filename : pars->filenames) {
        pars->imgs.push_back(yimage());
        auto& img = pars->imgs.back();
        img.filename = filename;
        if (yimg::is_hdr_filename(filename)) {
            img.hdr = yimg::load_image4f(filename);
        } else {
            img.ldr = yimg::load_image4b(filename);
        }
        if (!img.hdr && !img.ldr) {
            log_fatal("cannot load image %s\n", img.filename.c_str());
        }
        img.tex_glid = 0;
    }
    return true;
}

void parse_cmdline(params* pars, int argc, char** argv) {
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"none", (int)ym::tonemap_type::none},
        {"srgb", (int)ym::tonemap_type::srgb},
        {"gamma", (int)ym::tonemap_type::gamma},
        {"filmic", (int)ym::tonemap_type::filmic}};

    auto parser =
        yu::cmdline::make_parser(argc, argv, "yimview", "view images");
    pars->exposure =
        parse_optf(parser, "--exposure", "-e", "hdr image exposure", 0);
    pars->gamma = parse_optf(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    pars->tonemap = (ym::tonemap_type)parse_opte(parser, "--tonemap", "-t",
        "hdr image tonemap", (int)ym::tonemap_type::srgb, tmtype_names);
    pars->filenames = parse_argas(parser, "image", "image filename", {}, true);

    // done
    check_parser(parser);
}

void text_callback(ygui::window* win, unsigned int key) {
    auto pars = (params*)get_user_pointer(win);
    switch (key) {
        case ' ':
        case '.':
            pars->cur_img = (pars->cur_img + 1) % pars->imgs.size();
            break;
        case ',':
            pars->cur_img = (pars->cur_img - 1 + (int)pars->imgs.size()) %
                            pars->imgs.size();
            break;
        case '-':
        case '_': pars->zoom /= 2; break;
        case '+':
        case '=': pars->zoom *= 2; break;
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '1':
            pars->exposure = 0;
            pars->gamma = 1;
            break;
        case '2':
            pars->exposure = 0;
            pars->gamma = 2.2f;
            break;
        case 'z': pars->zoom = 1; break;
        case 'h':
            // TODO: hud
            break;
        default: log_fatal("unsupported key\n"); break;
    }
}

void draw_image(ygui::window* win) {
    auto pars = (params*)get_user_pointer(win);

    auto& img = pars->imgs[pars->cur_img];

    // begin frame
    yglu::clear_buffers(
        {pars->background, pars->background, pars->background, 0});

    // draw image
    auto ws = ygui::get_widget_size(win);
    auto wwh = ygui::get_window_size(win);
    auto fwh = ygui::get_framebuffer_size(win);
    fwh.x -= (ws * fwh.x) / wwh.x;
    wwh.x -= ws;
    yglu::set_viewport({0, 0, fwh.x, fwh.y});
    yglu::shade_image(img.tex_glid, wwh.x, wwh.y, pars->offset.x,
        pars->offset.y, pars->zoom, pars->tonemap, pars->exposure, pars->gamma);
}

template <typename T>
ym::vec<T, 4> lookup_image(
    int w, int h, int nc, const T* pixels, int x, int y, T one) {
    if (x < 0 || y < 0 || x > w - 1 || y > h - 1) return {0, 0, 0, 0};
    auto v = ym::vec<T, 4>{0, 0, 0, 0};
    auto vv = pixels + ((w * y) + x) * nc;
    switch (nc) {
        case 1: v = {vv[0], 0, 0, one}; break;
        case 2: v = {vv[0], vv[1], 0, one}; break;
        case 3: v = {vv[0], vv[1], vv[2], one}; break;
        case 4: v = {vv[0], vv[1], vv[2], vv[3]}; break;
        default: assert(false);
    }
    return v;
}

void draw_widgets(ygui::window* win) {
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"none", (int)ym::tonemap_type::none},
        {"srgb", (int)ym::tonemap_type::srgb},
        {"gamma", (int)ym::tonemap_type::gamma},
        {"filmic", (int)ym::tonemap_type::filmic}};

    auto pars = (params*)get_user_pointer(win);
    auto& img = pars->imgs[pars->cur_img];
    auto mouse_pos = (ym::vec2f)get_mouse_posf(win);
    if (begin_widgets(win, "yimview")) {
        label_widget(win, "filename", img.filename);
        label_widget(win, "w", img.width());
        label_widget(win, "h", img.height());
        if (img.hdr) {
            combo_widget(win, "tonemap", (int*)&pars->tonemap, tmtype_names);
            slider_widget(win, "exposure", &pars->exposure, -20, 20, 1);
            slider_widget(win, "gamma", &pars->gamma, 0.1, 5, 0.1);
        }
        auto offset = ym::vec2i{(int)pars->offset.x, (int)pars->offset.y};
        slider_widget(win, "offset", &offset, -100, 100, 1);
        pars->offset = ym::vec2f{(float)offset.x, (float)offset.y};
        slider_widget(win, "zoom", &pars->zoom, 0.01, 10, 0.1);
        auto xy = (mouse_pos - pars->offset) / pars->zoom;
        auto ij = ym::vec2i{(int)round(xy[0]), (int)round(xy[1])};
        auto inside = ij[0] >= 0 && ij[1] >= 0 && ij[0] < img.width() &&
                      ij[1] < img.height();
        if (img.hdr) {
            auto p = (inside) ? img.hdr[ij] : ym::zero4f;
            label_widget(win, "r", p.x);
            label_widget(win, "g", p.y);
            label_widget(win, "b", p.z);
            label_widget(win, "a", p.w);
        }
        if (img.ldr) {
            auto p = (inside) ? img.ldr[ij] : ym::zero4b;
            label_widget(win, "r", p.x);
            label_widget(win, "g", p.y);
            label_widget(win, "b", p.z);
            label_widget(win, "a", p.w);
        }
    }
    end_widgets(win);
}

void window_refresh_callback(ygui::window* win) {
    draw_image(win);
    draw_widgets(win);
    swap_buffers(win);
}

void run_ui(params* pars) {
    // window
    auto win = ygui::init_window(
        pars->imgs[0].width() + 320, pars->imgs[0].height(), "yimview", pars);
    set_callbacks(win, text_callback, nullptr, window_refresh_callback);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;

    init_widgets(win);

    // load textures
    for (auto& img : pars->imgs) {
        if (img.hdr) {
            img.tex_glid = yglu::make_texture(img.hdr, false, false, true);
        } else if (img.ldr) {
            img.tex_glid = yglu::make_texture(img.ldr, false, false, true);
        }
    }

    while (!should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = get_mouse_posf(win);
        mouse_button = get_mouse_button(win);

        auto& img = pars->imgs[pars->cur_img];
        set_window_title(win,
            ("yimview | " + img.filename + " | " + std::to_string(img.width()) +
                "x" + std::to_string(img.height()))
                .c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !get_widget_active(win)) {
            switch (mouse_button) {
                case 1: pars->offset += mouse_pos - mouse_last; break;
                case 2:
                    pars->zoom *=
                        powf(2, (mouse_pos[0] - mouse_last[0]) * 0.001f);
                    break;
                default: break;
            }
        }

        // draw
        draw_image(win);
        draw_widgets(win);

        // swap buffers
        swap_buffers(win);

        // event hadling
        wait_events(win);
    }

    clear_widgets(win);
    clear_window(win);
}

int main(int argc, char* argv[]) {
    // command line params
    auto pars = new params();
    parse_cmdline(pars, argc, argv);

    // loading images
    if (!load_images(pars)) return 1;

    // run ui
    run_ui(pars);

    // done
    delete pars;
    return 0;
}
