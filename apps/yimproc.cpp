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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"
#include "CLI11.hpp"
using namespace std::literals;

#if 0
template <typename Image>
Image make_image_grid(const std::vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width * tilex;
    auto height = imgs[0].height * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = ygl::make_image4f(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.size != imgs[0].size) {
            ygl::log_fatal("images of different sizes are not accepted");
        }
        auto ox = (img_idx % tilex) * img.width,
             oy = (img_idx / tilex) * img.height;
        if (ret.hdr) {
            for (auto j = 0; j < img.height; j++) {
                for (auto i = 0; i < img.width; i++) {
                    ret.hdr[{i + ox, j + oy}] = img.hdr.at(i, j);
                }
            }
        } else {
            for (auto j = 0; j < img.height; j++) {
                for (auto i = 0; i < img.width; i++) {
                    ret.ldr[{i + ox, j + oy}] = img.ldr.at(i, j);
                }
            }
        }
    }
    return ret;
}
#endif

ygl::image4f filter_bilateral(const ygl::image4f& img, float spatial_sigma,
    float range_sigma, const std::vector<ygl::image4f>& features,
    const std::vector<float>& features_sigma) {
    auto filtered = ygl::make_image4f(img.width, img.height);
    auto filter_width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = std::vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.height; j++) {
        for (auto i = 0; i < img.width; i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -filter_width; fj <= filter_width; fj++) {
                for (auto fi = -filter_width; fi <= filter_width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width || jj >= img.height) continue;
                    auto uv = ygl::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img.at(i, j) - img.at(ii, jj);
                    auto w = (float)std::exp(-dot(uv, uv) * sw) *
                             (float)std::exp(-dot(rgb, rgb) * rw);
                    for (auto fi = 0; fi < features.size(); fi++) {
                        auto feat =
                            features[fi].at(i, j) - features[fi].at(ii, jj);
                        w *= exp(-dot(feat, feat) * fw[fi]);
                    }
                    av += w * img.at(ii, jj);
                    aw += w;
                }
            }
            filtered.at(i, j) = av / aw;
        }
    }
    return filtered;
}

ygl::image4f filter_bilateral(
    const ygl::image4f& img, float spatial_sigma, float range_sigma) {
    auto filtered = ygl::make_image4f(img.width, img.height);
    auto fwidth = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < img.height; j++) {
        for (auto i = 0; i < img.width; i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -fwidth; fj <= fwidth; fj++) {
                for (auto fi = -fwidth; fi <= fwidth; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width || jj >= img.height) continue;
                    auto uv = ygl::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img.at(i, j) - img.at(ii, jj);
                    auto w = std::exp(-dot(uv, uv) * sw) *
                             std::exp(-dot(rgb, rgb) * rw);
                    av += w * img.at(ii, jj);
                    aw += w;
                }
            }
            filtered.at(i, j) = av / aw;
        }
    }
    return filtered;
}

int main(int argc, char* argv[]) {
    // command line parameters
    auto filename = "img.png"s;      // input image
    auto output = "out.png"s;        // output image
    auto tonemap = false;            // enable tonemapping
    auto exposure = 0.0f;            // tonemap exposure
    auto gamma = 1.0f;               // tonemap gamma
    auto filmic = false;             // tonemap filmic
    auto res_width = 0;              // resize width
    auto res_height = 0;             // resize height
    auto alpha_filename = ""s;       // file to copy alpha from
    auto coloralpha_filename = ""s;  // file to set alpha from color
    auto spatial_sigma = 0.0f;       // spatial sigma for bilateral blur
    auto range_sigma = 0.0f;         // range sigma for bilateral blur

    // command line params
    CLI::App parser("image processing utility", "yimproc");
    parser.add_flag("--tonemap,-t", tonemap, "Tonemap image");
    parser.add_option("--exposure,-e", exposure, "Tonemap exposure");
    parser.add_option("--gamma,-g", gamma, "Tonemap gamma.");
    parser.add_flag("--filmic,-f", filmic, "Tonemap uses filmic curve");
    parser.add_option(
        "--res-width", res_width, "resize width (0 to maintain aspect)");
    parser.add_option(
        "--res-height", res_height, "resize height (0 to maintain aspect)");
    parser.add_option("--spatial-sigma", spatial_sigma, "blur spatial sigma");
    parser.add_option(
        "--range-sigma", range_sigma, "bilateral blur range sigma");
    parser.add_option(
        "--set-alpha", alpha_filename, "set alpha as this image alpha");
    parser.add_option("--set-color-as-alpha", coloralpha_filename,
        "set alpha as this image color");
    parser.add_option("--output,-o", output, "output image filename");
    parser.add_option("filename", filename, "input image filename")
        ->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // load
    auto img = ygl::image4f();
    try {
        img = ygl::load_image(filename);
    } catch (std::exception& e) {
        printf("cannot load image %s\n", filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }

    // set alpha
    if (alpha_filename != "") {
        auto alpha = ygl::load_image(alpha_filename);
        if (img.width != alpha.width || img.height != alpha.height) {
            printf("bad image size\n");
            exit(1);
        }
        for (auto j = 0; j < img.height; j++)
            for (auto i = 0; i < img.width; i++)
                img.at(i, j).w = alpha.at(i, j).w;
    }

    // set alpha
    if (coloralpha_filename != "") {
        auto alpha = ygl::load_image(coloralpha_filename);
        if (img.width != alpha.width || img.height != alpha.height) {
            printf("bad image size\n");
            exit(1);
        }
        for (auto j = 0; j < img.height; j++)
            for (auto i = 0; i < img.width; i++)
                img.at(i, j).w = ygl::luminance(alpha.at(i, j));
    }

    // resize
    if (res_width || res_height) {
        img = resize_image(img, res_width, res_height);
    }

    // bilateral
    if (spatial_sigma && range_sigma) {
        img = filter_bilateral(img, spatial_sigma, range_sigma, {}, {});
    }

    // hdr correction
    if (tonemap) img = ygl::tonemap_image(img, exposure, gamma, filmic);

    // save
    try {
        ygl::save_image(output, img);
    } catch (std::exception& e) {
        printf("cannot save image %s\n", output.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }

    // done
    return 0;
}
