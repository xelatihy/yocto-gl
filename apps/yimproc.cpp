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

#define YCMD_INLINE
#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_math.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../yocto/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../yocto/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "../yocto/stb_image_resize.h"

using byte = unsigned char;

struct image {
    int width = 0, height = 0, ncomp = 0;
    float *hdr = nullptr;
    byte *ldr = nullptr;

    ~image() {
        if (ldr) delete ldr;
        if (hdr) delete hdr;
    }
};

std::shared_ptr<image> load_image(const std::string &filename) {
    auto img = std::make_shared<image>();
    if (ycmd::get_extension(filename) == ".hdr") {
        img->hdr = stbi_loadf(filename.c_str(), &img->width, &img->height,
                              &img->ncomp, 0);
    } else {
        img->ldr = stbi_load(filename.c_str(), &img->width, &img->height,
                             &img->ncomp, 0);
    }
    if (!img->ldr && !img->hdr)
        throw std::runtime_error("cannot load image " + filename);
    return img;
}

void save_image(const std::string &filename,
                const std::shared_ptr<image> &img) {
    if (ycmd::get_extension(filename) == ".hdr") {
        if (!img->hdr) throw std::invalid_argument("hdr data required");
        stbi_write_hdr(filename.c_str(), img->width, img->height, img->ncomp,
                       img->hdr);
    } else if (ycmd::get_extension(filename) == ".png") {
        if (!img->ldr) throw std::invalid_argument("ldr data required");
        stbi_write_png(filename.c_str(), img->width, img->height, img->ncomp,
                       img->ldr, img->width * img->ncomp);
    } else {
        throw std::invalid_argument("unsupported output extension " + filename);
    }
}

std::shared_ptr<image> resize_image(const std::shared_ptr<image> &img,
                                    int width, int height) {
    if (width < 0 && height < 0)
        throw std::invalid_argument("at least argument should be >0");
    auto res = std::make_shared<image>();
    res->width =
        (width > 0)
            ? width
            : (int)std::round(img->width * (height / (float)img->height));
    res->height =
        (height < 0)
            ? (int)std::round(img->height * (width / (float)img->width))
            : height;
    res->ncomp = img->ncomp;
    if (img->hdr) {
        res->hdr = new float[img->width * img->height * img->ncomp];
        auto img_stride = sizeof(float) * img->width * img->ncomp;
        auto res_stride = sizeof(float) * res->width * res->ncomp;
        stbir_resize_float(img->hdr, img->width, img->height, img_stride,
                           res->hdr, res->width, res->height, res_stride,
                           img->ncomp);
    } else {
        res->ldr = new byte[img->width * img->height * img->ncomp];
        auto img_stride = sizeof(byte) * img->width * img->ncomp;
        auto res_stride = sizeof(byte) * res->width * res->ncomp;
        stbir_resize_uint8_srgb(
            img->ldr, img->width, img->height, img_stride, res->ldr, res->width,
            res->height, res_stride, img->ncomp,
            (img->ncomp == 4) ? 3 : STBIR_ALPHA_CHANNEL_NONE, 0);
    }
    return res;
}

std::shared_ptr<image> make_image_grid(
    const std::vector<std::shared_ptr<image>> &imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto ret = std::make_shared<image>();
    ret->width = imgs[0]->width * tilex;
    ret->height = imgs[0]->height * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    ret->ncomp = imgs[0]->ncomp;
    auto nvalues = ret->width * ret->height * ret->ncomp;
    if (imgs[0]->hdr) {
        ret->hdr = new float[nvalues];
        memset(ret->hdr, 0, sizeof(float) * nvalues);
    } else {
        ret->ldr = new byte[nvalues];
        memset(ret->ldr, 0, sizeof(byte) * nvalues);
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

std::shared_ptr<image> make_image_grid(
    const std::vector<std::shared_ptr<image>> &imgs, int tilex, int width,
    int height) {
    auto resized = std::vector<std::shared_ptr<image>>();
    for (auto img : imgs) resized.push_back(resize_image(img, width, height));
    return make_image_grid(resized, tilex);
}

int main(int argc, char *argv[]) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "process images");
    auto command = ycmd::parse_args(parser, "command", "command to execute", "",
                                    true, {"resize"});
    auto output = ycmd::parse_opts(parser, "--output", "-o",
                                   "output image filename", "", true);
    auto width = ycmd::parse_opti(parser, "--width", "",
                                  "width (-1 to maintain aspect)", -1);
    auto height = ycmd::parse_opti(parser, "--height", "",
                                   "height (-1 to maintain aspect)", -1);
    auto filenames = ycmd::parse_argas(
        parser, "filenames", "input image filenames", {}, true, -1, {});

    // done command line parameters
    ycmd::check_parser(parser);

    // load images
    std::vector<std::shared_ptr<image>> imgs;
    for (auto filename : filenames) imgs.push_back(load_image(filename));

    // decalre output
    auto out = std::shared_ptr<image>();

    // switch on commands
    if (command == "resize") {
        out = resize_image(imgs[0], width, height);
    }

    // save output
    save_image(output, out);

    // done
    return 0;
}
