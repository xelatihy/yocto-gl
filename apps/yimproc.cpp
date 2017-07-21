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

#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_utils.h"

#include "string.h"

// typedef
using byte = unsigned char;

//
// Simple image structure to ease the passing of parameters.
//
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
};

//
// Initializes an image
//
yimage make_image(int width, int height, bool hdr) {
    auto img = yimage();
    if (hdr)
        img.hdr = ym::image<ym::vec4f>(width, height);
    else
        img.ldr = ym::image<ym::vec4b>(width, height);
    return img;
}

//
// Loads an image
//
yimage load_image(const std::string& filename) {
    auto img = yimage();
    auto ext = yu::path::get_extension(filename);
    if (ext == ".hdr") {
        img.hdr = yimg::load_image4f(filename);
    } else {
        img.ldr = yimg::load_image4b(filename);
    }
    if (!img.ldr && !img.hdr) {
        printf("could not load image %s\f", filename.c_str());
        exit(0);
    }
    return img;
}

//
// Saves an image
//
void save_image(const std::string& filename, const yimage& img) {
    if (img.ldr) { yimg::save_image4b(filename, img.ldr); }
    if (img.hdr) { yimg::save_image4f(filename, img.hdr); }
}

//
// Tone mapping HDR to LDR images.
//
yimage tonemap_image(
    const yimage& hdr, float exposure, ym::tonemap_type tm, float gamma) {
    if (!hdr.hdr) {
        printf("tonemap hdr only\n");
        exit(1);
    }
    auto ldr = make_image(hdr.width(), hdr.height(), false);
    tonemap_image(hdr.hdr, ldr.ldr, tm, exposure, gamma);
    return ldr;
}

//
// Resize image.
//
yimage resize_image(const yimage& img, int res_width, int res_height) {
    if (res_width < 0 && res_height < 0)
        throw std::invalid_argument("at least argument should be >0");
    if (res_width < 0)
        res_width =
            (int)std::round(img.width() * (res_height / (float)img.height()));
    if (res_height < 0)
        res_height =
            (int)std::round(img.height() * (res_width / (float)img.width()));
    auto res = make_image(res_width, res_height, (bool)img.hdr);
    if (img.hdr) yimg::resize_image(img.hdr, res.hdr);
    if (img.ldr) yimg::resize_image(img.ldr, res.ldr);
    return res;
}

yimage make_image_grid(const std::vector<yimage>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width() * tilex;
    auto height =
        imgs[0].height() * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = make_image(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.width() != imgs[0].width() ||
            img.height() != imgs[0].height()) {
            printf("images of different sizes are not accepted\n");
            exit(1);
        }
        auto ox = (img_idx % tilex) * img.width(),
             oy = (img_idx / tilex) * img.height();
        if (ret.hdr) {
            for (auto j = 0; j < img.height(); j++) {
                for (auto i = 0; i < img.width(); i++) {
                    ret.hdr[{i + ox, j + oy}] = img.hdr[{i, j}];
                }
            }
        } else {
            for (auto j = 0; j < img.height(); j++) {
                for (auto i = 0; i < img.width(); i++) {
                    ret.ldr[{i + ox, j + oy}] = img.ldr[{i, j}];
                }
            }
        }
    }
    return ret;
}

yimage make_image_grid(
    const std::vector<yimage>& imgs, int tilex, int width, int height) {
    auto resized = std::vector<yimage>();
    for (auto img : imgs) resized.push_back(resize_image(img, width, height));
    return make_image_grid(resized, tilex);
}

int main(int argc, char* argv[]) {
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"none", (int)ym::tonemap_type::none},
        {"srgb", (int)ym::tonemap_type::srgb},
        {"gamma", (int)ym::tonemap_type::gamma},
        {"filmic", (int)ym::tonemap_type::filmic}};

    // command line params
    auto parser = yu::cmdline::make_parser(argc, argv, "process images");
    auto command = parse_args(parser, "command", "command to execute", "", true,
        {"resize", "tonemap"});
    auto output =
        parse_opts(parser, "--output", "-o", "output image filename", "", true);
    auto width = parse_opti(
        parser, "--width", "-w", "width (-1 to maintain aspect)", -1);
    auto height = parse_opti(
        parser, "--height", "-h", "height (-1 to maintain aspect)", -1);
    auto exposure = parse_optf(parser, "--exposure", "-e", "hdr exposure", 0);
    auto gamma = parse_optf(parser, "--gamma", "-g", "hdr gamma", 2.2f);
    auto tonemap = (ym::tonemap_type)parse_opte(parser, "--tonemap", "-t",
        "hdr tonemap", (int)ym::tonemap_type::srgb, tmtype_names);
    auto filenames = parse_argas(
        parser, "filenames", "input image filenames", {}, true, -1, {});

    // done command line parameters
    check_parser(parser);

    // load images
    std::vector<yimage> imgs;
    for (auto filename : filenames) imgs.push_back(load_image(filename));

    // decalre output
    auto out = yimage();

    // switch on commands
    if (command == "resize") { out = resize_image(imgs[0], width, height); }
    if (command == "tonemap") {
        out = tonemap_image(imgs[0], exposure, tonemap, gamma);
    }

    // save output
    save_image(output, out);

    // done
    return 0;
}
