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
// Check if hdr
//
bool is_hdr(const std::string& filename) {
    return yimg::is_hdr_filename(filename);
}

//
// exit with an error
//
void exit_error(const std::string& err) {
    printf("error: %s\n", err.c_str());
    exit(1);
}

//
// Load hdr
//
ym::image4f load_hdr(const std::string& filename) {
    auto img = yimg::load_image4f(filename);
    if (!img) exit_error("cannot load image " + filename);
    return img;
}

//
// Load ldr
//
ym::image4b load_ldr(const std::string& filename) {
    auto img = yimg::load_image4b(filename);
    if (!img) exit_error("cannot load image " + filename);
    return img;
}

//
// Save hdr
//
void save_hdr(const std::string& filename, const ym::image4f& img) {
    if (!yimg::save_image4f(filename, img))
        exit_error("cannot save image " + filename);
}

//
// Save ldr
//
void save_ldr(const std::string& filename, const ym::image4b& img) {
    if (!yimg::save_image4b(filename, img))
        exit_error("cannot save image " + filename);
}

//
// Resize image.
//
template <typename T>
ym::image<T> resize_image(
    const ym::image<T>& img, int res_width, int res_height) {
    if (res_width < 0 && res_height < 0)
        exit_error("at least argument should be >0");
    if (res_width < 0)
        res_width =
            (int)std::round(img.width() * (res_height / (float)img.height()));
    if (res_height < 0)
        res_height =
            (int)std::round(img.height() * (res_width / (float)img.width()));
    auto res = ym::image<T>(res_width, res_height);
    yimg::resize_image(img, res);
    return res;
}

template <typename T>
ym::image<T> make_image_grid(const std::vector<ym::image<T>>& imgs, int tilex) {
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

#if 0
yimage make_image_grid(
    const std::vector<yimage>& imgs, int tilex, int width, int height) {
    auto resized = std::vector<yimage>();
    for (auto img : imgs) resized.push_back(resize_image(img, width, height));
    return make_image_grid(resized, tilex);
}
#endif

ym::image4f filter_bilateral(const ym::image4f& img, float spatial_sigma,
    float range_sigma, const std::vector<ym::image4f>& features,
    const std::vector<float>& features_sigma) {
    auto filtered = ym::image4f(img.width(), img.height());
    auto width = (int)ym::ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = std::vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = ym::zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
                    auto uv = ym::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[{i, j}] - img[{ii, jj}];
                    auto w = ym::exp(-lengthsqr(uv) * sw) *
                             ym::exp(-lengthsqr(rgb) * rw);
                    for (auto fi = 0; fi < features.size(); fi++) {
                        auto feat =
                            features[fi][{i, j}] - features[fi][{ii, jj}];
                        w *= ym::exp(-lengthsqr(feat) * fw[fi]);
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

ym::image4f filter_bilateral(
    const ym::image4f& img, float spatial_sigma, float range_sigma) {
    auto filtered = ym::image4f(img.width(), img.height());
    auto width = (int)ym::ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = ym::zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
                    auto uv = ym::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[{i, j}] - img[{ii, jj}];
                    auto w = ym::exp(-lengthsqr(uv) * sw) *
                             ym::exp(-lengthsqr(rgb) * rw);
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
    static auto tmtype_names =
        std::vector<std::pair<std::string, ym::tonemap_type>>{
            {"none", ym::tonemap_type::none}, {"srgb", ym::tonemap_type::srgb},
            {"gamma", ym::tonemap_type::gamma},
            {"filmic", ym::tonemap_type::filmic}};

    // command line params
    auto parser = yu::cmdline::make_parser(argc, argv, "process images");
    auto command = parse_args(parser, "command", "command to execute", "", true,
        {"resize", "tonemap", "bilateral"});
    auto output =
        parse_opts(parser, "--output", "-o", "output image filename", "", true);
    if (command == "resize") {
        auto width = parse_opti(
            parser, "--width", "-w", "width (-1 to maintain aspect)", -1);
        auto height = parse_opti(
            parser, "--height", "-h", "height (-1 to maintain aspect)", -1);
        auto filename =
            parse_args(parser, "filename", "input image filename", "", true);
        check_parser(parser);

        if (is_hdr(filename)) {
            auto img = load_hdr(filename);
            auto out = resize_image(img, width, height);
            save_hdr(output, out);
        } else {
            auto img = load_ldr(filename);
            auto out = resize_image(img, width, height);
            save_ldr(output, out);
        }
    } else if (command == "tonemap") {
        auto exposure =
            parse_optf(parser, "--exposure", "-e", "hdr exposure", 0);
        auto gamma = parse_optf(parser, "--gamma", "-g", "hdr gamma", 2.2f);
        auto tonemap = parse_opte(parser, "--tonemap", "-t", "hdr tonemap",
            ym::tonemap_type::srgb, tmtype_names);
        auto filename =
            parse_args(parser, "filename", "input image filename", "", true);
        check_parser(parser);

        auto img = load_hdr(filename);
        auto out = ym::tonemap_image(img, tonemap, exposure, gamma);
        save_ldr(output, out);
    } else if (command == "bilateral") {
        auto spatial_sigma =
            parse_opti(parser, "--spatial-sigma", "-s", "spatial sigma", 3);
        auto range_sigma =
            parse_optf(parser, "--range-sigma", "-r", "range sigma", 0.1f);
        auto feature_sigma = parse_optf(
            parser, "--features-sigma", "-f", "features sigmas", 0.1f);
        auto filename =
            parse_args(parser, "filename", "input image filename", "", true);
        auto ffilenames = parse_argas(
            parser, "features", "input features filename", {}, -1, false);
        check_parser(parser);

        auto img = load_hdr(filename);
        auto features = std::vector<ym::image4f>();
        auto features_sigma = std::vector<float>();
        for (auto ffilename : ffilenames) {
            features.push_back(load_hdr(ffilename));
            features_sigma.push_back(feature_sigma);
        }
        auto out = filter_bilateral(
            img, spatial_sigma, range_sigma, features, features_sigma);
        save_hdr(output, out);
    } else {
        check_parser(parser);
    }

    // done
    return 0;
}
