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

struct app_state {
    std::string filename = "img.png";
    std::string output = "out.png";
    float exposure = 0.0f;
    float gamma = 2.2f;
    bool filmic = false;
    int resize_width = 0;
    int resize_height = 0;
    ygl::vec4f multiply_color = {1, 1, 1, 1};
    std::string alpha_filename = ""s;
    std::string coloralpha_filename = ""s;

    float spatial_sigma = 0.0f;
    float range_sigma = 0.0f;

    ygl::image4f img = {};
};

#if 0
template <typename Image>
Image make_image_grid(const std::vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width * tilex;
    auto height = imgs[0].height * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = ygl::make_image4f(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.width() != imgs[0].width || img.height() != imgs[0].height) {
            ygl::log_fatal("images of different sizes are not accepted");
        }
        auto ox = (img_idx % tilex) * img.width(),
             oy = (img_idx / tilex) * img.height();
        if (ret.hdr) {
            for (auto j = 0; j < img.height(); j++) {
                for (auto i = 0; i < img.width(); i++) {
                    ret.hdr[{i + ox, j + oy}] = img.hdr.at(i, j);
                }
            }
        } else {
            for (auto j = 0; j < img.height(); j++) {
                for (auto i = 0; i < img.width(); i++) {
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
    auto filtered = ygl::image4f{img.width(), img.height()};
    auto filter_width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = std::vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -filter_width; fj <= filter_width; fj++) {
                for (auto fi = -filter_width; fi <= filter_width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
                    auto uv = ygl::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[{i, j}] - img[{ii, jj}];
                    auto w = (float)std::exp(-dot(uv, uv) * sw) *
                             (float)std::exp(-dot(rgb, rgb) * rw);
                    for (auto fi = 0; fi < features.size(); fi++) {
                        auto feat =
                            features[fi][{i, j}] - features[fi][{ii, jj}];
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

std::vector<ygl::vec4f> filter_bilateral(int width, int height,
    const std::vector<ygl::vec4f>& img, float spatial_sigma,
    float range_sigma) {
    auto filtered = std::vector<ygl::vec4f>(width * height);
    auto fwidth = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -fwidth; fj <= fwidth; fj++) {
                for (auto fi = -fwidth; fi <= fwidth; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= width || jj >= height) continue;
                    auto uv = ygl::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[i + width * j] - img[ii + width * jj];
                    auto w = std::exp(-dot(uv, uv) * sw) *
                             std::exp(-dot(rgb, rgb) * rw);
                    av += w * img[ii + width * jj];
                    aw += w;
                }
            }
            filtered[i + width * j] = av / aw;
        }
    }
    return filtered;
}

int main(int argc, char* argv[]) {
    auto app = std::make_shared<app_state>();

    // command line params
    CLI::App parser("image processing utility", "yimproc");
    parser.add_option("--exposure,-e", app->exposure, "Hdr exposure");
    parser.add_option("--gamma,-g", app->gamma, "Display gamma.");
    parser.add_flag("--filmic,-f", app->filmic, "filmic tone mapping");
    parser.add_option("--resize-width,-w", app->resize_width,
        "width (0 to maintain aspect)", 0);
    parser.add_option("--resize-height,-h", app->resize_height,
        "height (0 to maintain aspect)", 0);
    parser.add_option(
        "--multiply-color", app->multiply_color, "multiply by this color");
    parser.add_option(
        "--spatial-sigma,-s", app->spatial_sigma, "blur spatial sigma");
    parser.add_option(
        "--range-sigma,-r", app->range_sigma, "bilateral blur range sigma");
    parser.add_option(
        "--set-alpha", app->alpha_filename, "set alpha as this image alpha");
    parser.add_option("--set-color-as-alpha", app->coloralpha_filename,
        "set alpha as this image color");
    parser.add_option("--output,-o", app->output, "output image filename");
    parser.add_option("filename", app->filename, "input image filename")
        ->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // load
    try {
        app->img = ygl::load_image(app->filename, app->gamma);
    } catch (std::exception& e) {
        std::cout << "cannot load image" << app->filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // set alpha
    if (app->alpha_filename != "") {
        auto alpha = ygl::load_image(app->alpha_filename, app->gamma);
        if (app->img.width() != alpha.width() ||
            app->img.height() != alpha.height()) {
            std::cout << "bad image size\n";
            exit(1);
        }
        for (auto i = 0; i < app->img.size(); i++) app->img[i].w = alpha[i].w;
    }

    // set alpha
    if (app->coloralpha_filename != "") {
        auto alpha = ygl::load_image(app->coloralpha_filename, app->gamma);
        if (app->img.width() != alpha.width() ||
            app->img.height() != alpha.height()) {
            std::cout << "bad image size\n";
            exit(1);
        }
        for (auto i = 0; i < app->img.size(); i++) {
            auto& p = alpha[i];
            app->img[i].w = (p.x + p.y + p.z) / 3;
        }
    }

    // multiply
    if (app->multiply_color != ygl::vec4f{1, 1, 1, 1}) {
        for (auto& c : app->img) c *= app->multiply_color;
    }

    // resize
    if (app->resize_width || app->resize_height) {
        app->img =
            resize_image(app->img, app->resize_width, app->resize_height);
    }

    // bilateral
    if (app->spatial_sigma && app->range_sigma) {
        app->img = filter_bilateral(
            app->img, app->spatial_sigma, app->range_sigma, {}, {});
    }

    // exposure
    if (app->exposure) app->img = expose_image(app->img, app->exposure);

    // filmic tone transformations
    if (app->filmic) app->img = filmic_tonemap_image(app->img);

    // save
    try {
        ygl::save_image(app->output, app->img, app->gamma);
    } catch (std::exception& e) {
        std::cout << "cannot save image" << app->output << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // done
    return 0;
}
