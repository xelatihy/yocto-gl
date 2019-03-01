//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

#if 0
template <typename Image>
Image make_image_grid(const vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].size().x * tilex;
    auto height = imgs[0].size().y * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = init_image(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (extents(img) != extents(imgs[0])) {
            exit_error("images of different sizes are not accepted");
        }
        auto ox = (img_idx % tilex) * img.size().x,
             oy = (img_idx / tilex) * img.size().y;
        if (ret.hdr) {
            for (auto j = 0; j < img.size().y; j++) {
                for (auto i = 0; i < img.size().x; i++) {
                    ret.hdr[{i + ox, j + oy}] = img.hdr[{i,j}];
                }
            }
        } else {
            for (auto j = 0; j < img.size().y; j++) {
                for (auto i = 0; i < img.size().x; i++) {
                    ret.ldr[{i + ox, j + oy}] = img.ldr[{i,j}];
                }
            }
        }
    }
    return ret;
}
#endif

image4f filter_bilateral(const image4f& img, float spatial_sigma,
    float range_sigma, const vector<image4f>& features,
    const vector<float>& features_sigma) {
    auto filtered     = image{img.size(), zero4f};
    auto filter_width = (int)ceil(2.57f * spatial_sigma);
    auto sw           = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw           = 1 / (2.0f * range_sigma * range_sigma);
    auto fw           = vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto av = zero4f;
            auto aw = 0.0f;
            for (auto fj = -filter_width; fj <= filter_width; fj++) {
                for (auto fi = -filter_width; fi <= filter_width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.size().x || jj >= img.size().y) continue;
                    auto uv  = vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[{i, j}] - img[{i, j}];
                    auto w   = (float)exp(-dot(uv, uv) * sw) *
                             (float)exp(-dot(rgb, rgb) * rw);
                    for (auto fi = 0; fi < features.size(); fi++) {
                        auto feat = features[fi][{i, j}] - features[fi][{i, j}];
                        w *= exp(-dot(feat, feat) * fw[fi]);
                    }
                    av += w * img[{ii, jj}];
                    aw += w;
                }
            }
            filtered[{i, j}] = av / aw;
        }
    }
    return filtered;
}

image4f filter_bilateral(
    const image4f& img, float spatial_sigma, float range_sigma) {
    auto filtered = image{img.size(), zero4f};
    auto fwidth   = (int)ceil(2.57f * spatial_sigma);
    auto sw       = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw       = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto av = zero4f;
            auto aw = 0.0f;
            for (auto fj = -fwidth; fj <= fwidth; fj++) {
                for (auto fi = -fwidth; fi <= fwidth; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.size().x || jj >= img.size().y) continue;
                    auto uv  = vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[{i, j}] - img[{ii, jj}];
                    auto w = exp(-dot(uv, uv) * sw) * exp(-dot(rgb, rgb) * rw);
                    av += w * img[{ii, jj}];
                    aw += w;
                }
            }
            filtered[{i, j}] = av / aw;
        }
    }
    return filtered;
}

int main(int argc, char* argv[]) {
    // command line parameters
    auto tonemap             = false;
    auto exposure            = 0.0f;
    auto srgb                = true;
    auto filmic              = false;
    auto resize_width        = 0;
    auto resize_height       = 0;
    auto spatial_sigma       = 0.0f;
    auto range_sigma         = 0.0f;
    auto alpha_filename      = ""s;
    auto coloralpha_filename = ""s;
    auto output              = "out.png"s;
    auto filename            = "img.hdr"s;

    // parse command line
    auto parser = CLI::App{"Transform images"};
    parser.add_flag("--tonemap,!--no-tonemap,-t", tonemap, "Tonemap image");
    parser.add_option("--exposure,-e", exposure, "Tonemap exposure");
    parser.add_flag("--srgb,!--no-srgb", srgb, "Tonemap to sRGB.");
    parser.add_flag(
        "--filmic,!--no-filmic,-f", filmic, "Tonemap uses filmic curve");
    parser.add_option(
        "--resize-width", resize_width, "resize size (0 to maintain aspect)");
    parser.add_option(
        "--resize-height", resize_height, "resize size (0 to maintain aspect)");
    parser.add_option("--spatial-sigma", spatial_sigma, "blur spatial sigma");
    parser.add_option(
        "--range-sigma", range_sigma, "bilateral blur range sigma");
    parser.add_option(
        "--set-alpha", alpha_filename, "set alpha as this image alpha");
    parser.add_option("--set-color-as-alpha", coloralpha_filename,
        "set alpha as this image color");
    parser.add_option("--output,-o", output, "output image filename")
        ->required(true);
    parser.add_option("filename", filename, "input image filename")
        ->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // load
    auto img = image4f();
    try {
        load_image(filename, img);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // set alpha
    if (alpha_filename != "") {
        auto alpha = image4f();
        try {
            load_image(alpha_filename, alpha);
        } catch (const std::exception& e) {
            exit_error(e.what());
        }
        if (img.size() != alpha.size()) {
            exit_error("bad image size");
            exit(1);
        }
        for (auto j = 0; j < img.size().y; j++)
            for (auto i = 0; i < img.size().x; i++)
                img[{i, j}].w = alpha[{i, j}].w;
    }

    // set alpha
    if (coloralpha_filename != "") {
        auto alpha = image4f();
        try {
            load_image(coloralpha_filename, alpha);
        } catch (const std::exception& e) {
            exit_error(e.what());
        }
        if (img.size() != alpha.size()) {
            exit_error("bad image size");
            exit(1);
        }
        for (auto j = 0; j < img.size().y; j++)
            for (auto i = 0; i < img.size().x; i++)
                img[{i, j}].w = mean(xyz(alpha[{i, j}]));
    }

    // resize
    if (resize_width != 0 || resize_height != 0) {
        auto res = image4f{};
        resize_image(res, img, {resize_width, resize_height});
        img = res;
    }

    // bilateral
    if (spatial_sigma && range_sigma) {
        img = filter_bilateral(img, spatial_sigma, range_sigma, {}, {});
    }

    // hdr correction
    if (tonemap) {
        auto ldr = img;
        tonemap_image(img, ldr, exposure, filmic, srgb);
        img = ldr;
    }

    // save
    try {
        save_image(output, img);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // done
    return 0;
}
