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
struct simage {
    /// image width
    int width = 0;

    /// image height
    int height = 0;

    /// number of components
    int ncomp = 0;

    /// float data for hdr images if loaded
    float* hdr = nullptr;

    /// char data for ldr images if loaded
    byte* ldr = nullptr;

    /// default constructor
    simage() {}

    /// allocating constructor
    simage(int width, int height, int ncomp, bool ishdr)
        : width(width)
        , height(height)
        , ncomp(ncomp)
        , hdr((ishdr) ? new float[width * height * ncomp] : nullptr)
        , ldr((ishdr) ? nullptr : new byte[width * height * ncomp]) {}

    /// destructor
    ~simage() {
        if (ldr) delete[] ldr;
        if (hdr) delete[] hdr;
    }
};

//
// Initializes an image
//
simage* make_image(int width, int height, int ncomp, bool hdr) {
    return new simage(width, height, ncomp, hdr);
}

//
// Loads an image
//
simage* load_image(const std::string& filename) {
    auto img = new simage();
    auto ext = yu::path::get_extension(filename);
    if (ext == ".hdr") {
        img->hdr =
            yimg::load_imagef(filename, img->width, img->height, img->ncomp);

    } else {
        img->ldr =
            yimg::load_image(filename, img->width, img->height, img->ncomp);
    }
    return img;
}

//
// Saves an image
//
void save_image(const std::string& filename, const simage* img) {
    if (img->ldr) {
        yimg::save_image(
            filename, img->width, img->height, img->ncomp, img->ldr);
    } else if (img->hdr) {
        yimg::save_imagef(
            filename, img->width, img->height, img->ncomp, img->hdr);

    } else {
        assert(false);
    }
}

//
// Resize image.
//
simage* resize_image(const simage* img, int res_width, int res_height) {
    auto res = new simage();
    if (img->hdr) {
        yimg::resize_image(img->width, img->height, img->ncomp, img->hdr,
            res->width, res->height, res->hdr);
    } else if (img->ldr) {
        yimg::resize_image(img->width, img->height, img->ncomp, img->ldr,
            res->width, res->height, res->ldr);

    } else {
        assert(false);
    }
    return res;
}

//
// Tone mapping HDR to LDR images.
//
simage* tonemap_image(
    simage* hdr, float exposure, ym::tonemap_type tm, float gamma) {
    if (!hdr->hdr) throw std::invalid_argument("tonemap hdr only");
    auto ldr = make_image(hdr->width, hdr->height, hdr->ncomp, false);
    if (hdr->ncomp == 3) {
        tonemap_image(hdr->width, hdr->height, (ym::vec3f*)hdr->hdr,
            (ym::vec3b*)ldr->ldr, tm, exposure, gamma);
    } else {
        tonemap_image(hdr->width, hdr->height, (ym::vec4f*)hdr->hdr,
            (ym::vec4b*)ldr->ldr, tm, exposure, gamma);
    }
    return ldr;
}

//
// Resize image.
//
simage* resize_image(simage* img, int res_width, int res_height) {
    if (res_width < 0 && res_height < 0)
        throw std::invalid_argument("at least argument should be >0");
    if (res_width < 0)
        res_width =
            (int)std::round(img->width * (res_height / (float)img->height));
    if (res_height < 0)
        res_height =
            (int)std::round(img->height * (res_width / (float)img->width));
    auto res = make_image(res_width, res_height, img->ncomp, (bool)img->hdr);
    if (img->hdr)
        yimg::resize_image(img->width, img->height, img->ncomp, img->hdr,
            res_width, res_height, res->hdr);
    if (img->ldr)
        yimg::resize_image(img->width, img->height, img->ncomp, img->ldr,
            res_width, res_height, res->ldr);
    return res;
}

simage* make_image_grid(const std::vector<simage*>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto ret = new simage();
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

simage* make_image_grid(
    const std::vector<simage*>& imgs, int tilex, int width, int height) {
    auto resized = std::vector<simage*>();
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
    std::vector<simage*> imgs;
    for (auto filename : filenames) imgs.push_back(load_image(filename));

    // decalre output
    auto out = (simage*)nullptr;

    // switch on commands
    if (command == "resize") { out = resize_image(imgs[0], width, height); }
    if (command == "tonemap") {
        out = tonemap_image(imgs[0], exposure, tonemap, gamma);
    }

    // save output
    save_image(output, out);

    // cleanup
    for (auto img : imgs) delete img;
    delete out;

    // done
    return 0;
}
