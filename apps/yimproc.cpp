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

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"

yimg::simage *make_image_grid(const std::vector<yimg::simage *> &imgs,
                              int tilex) {
    auto nimgs = (int)imgs.size();
    auto ret = new yimg::simage();
    ret->width = imgs[0]->width * tilex;
    ret->height = imgs[0]->height * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    ret->ncomp = imgs[0]->ncomp;
    auto nvalues = ret->width * ret->height * ret->ncomp;
    if (imgs[0]->hdr) {
        ret->hdr = new float[nvalues];
        memset(ret->hdr, 0, sizeof(float) * nvalues);
    } else {
        ret->ldr = new yimg::byte[nvalues];
        memset(ret->ldr, 0, sizeof(yimg::byte) * nvalues);
    }
    auto img_idx = 0;
    for (auto img : imgs) {
        if (img->width != imgs[0]->width || img->height != imgs[0]->height ||
            img->ncomp != imgs[0]->ncomp) {
            throw std::invalid_argument(
                "images of different sizes are not accepted");
        }
        auto ox = (img_idx % tilex) * img->width,
             oy = (img_idx / tilex) * img->height;
        if (ret->hdr) {
            for (auto j = 0; j < img->height; j++) {
                for (auto i = 0; i < img->width; i++) {
                    auto retp = ret->hdr +
                                ((j + oy) * ret->width + i + ox) * ret->ncomp;
                    auto imgp = img->hdr + (j * img->width + i) * img->ncomp;
                    for (auto c = 0; c < ret->ncomp; c++) retp[c] = imgp[c];
                }
            }
        } else {
            for (auto j = 0; j < img->height; j++) {
                for (auto i = 0; i < img->width; i++) {
                    auto retp = ret->ldr +
                                ((j + oy) * ret->width + i + ox) * ret->ncomp;
                    auto imgp = img->ldr + (j * img->width + i) * img->ncomp;
                    for (auto c = 0; c < ret->ncomp; c++) retp[c] = imgp[c];
                }
            }
        }
    }
    return ret;
}

yimg::simage *make_image_grid(const std::vector<yimg::simage *> &imgs,
                              int tilex, int width, int height) {
    auto resized = std::vector<yimg::simage *>();
    for (auto img : imgs) resized.push_back(resize_image(img, width, height));
    return make_image_grid(resized, tilex);
}

int main(int argc, char *argv[]) {
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"default", (int)yimg::tonemap_type::def},
        {"linear", (int)yimg::tonemap_type::linear},
        {"srgb", (int)yimg::tonemap_type::srgb},
        {"gamma", (int)yimg::tonemap_type::gamma},
        {"filmic", (int)yimg::tonemap_type::filmic}};

    // command line params
    auto parser = ycmd::make_parser(argc, argv, "process images");
    auto command = ycmd::parse_args(parser, "command", "command to execute", "",
                                    true, {"resize", "tonemap"});
    auto output = ycmd::parse_opts(parser, "--output", "-o",
                                   "output image filename", "", true);
    auto width = ycmd::parse_opti(parser, "--width", "-w",
                                  "width (-1 to maintain aspect)", -1);
    auto height = ycmd::parse_opti(parser, "--height", "-h",
                                   "height (-1 to maintain aspect)", -1);
    auto exposure =
        ycmd::parse_optf(parser, "--exposure", "-e", "hdr exposure", 0);
    auto gamma = ycmd::parse_optf(parser, "--gamma", "-g", "hdr gamma", 2.2f);
    auto tonemap = (yimg::tonemap_type)ycmd::parse_opte(
        parser, "--tonemap", "-t", "hdr tonemap", (int)yimg::tonemap_type::srgb,
        tmtype_names);
    auto filenames = ycmd::parse_argas(
        parser, "filenames", "input image filenames", {}, true, -1, {});

    // done command line parameters
    ycmd::check_parser(parser);

    // load images
    std::vector<yimg::simage *> imgs;
    for (auto filename : filenames) imgs.push_back(yimg::load_image(filename));

    // decalre output
    auto out = (yimg::simage *)nullptr;

    // switch on commands
    if (command == "resize") {
        out = resize_image(imgs[0], width, height);
    }
    if (command == "tonemap") {
        out = yimg::tonemap_image(imgs[0], exposure, tonemap, gamma);
    }

    // save output
    save_image(output, out);

    // cleanup
    for (auto img : imgs) delete img;
    delete out;

    // done
    return 0;
}
