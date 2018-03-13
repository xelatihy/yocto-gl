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

#include "../yocto/yocto_gl.h"
using namespace std::literals;

// Load hdr
ygl::image4f load_hdr(const std::string& filename) {
    auto img = ygl::load_image4f(filename);
    if (img.empty()) ygl::log_fatal("cannot load image {}", filename);
    return img;
}

// Load ldr
ygl::image4b load_ldr(const std::string& filename) {
    auto img = ygl::load_image4b(filename);
    if (img.empty()) ygl::log_fatal("cannot load image {}", filename);
    return img;
}

// Save hdr
void save_hdr(const std::string& filename, const ygl::image4f& img) {
    if (!save_image4f(filename, img))
        ygl::log_fatal("cannot save image {}", filename);
}

// Save ldr
void save_ldr(const std::string& filename, const ygl::image4b& img) {
    if (!save_image4b(filename, img))
        ygl::log_fatal("cannot save image {}", filename);
}

// Resize image.
template <typename Image>
Image resize_image(const Image& img, int res_width, int res_height) {
    if (res_width < 0 && res_height < 0)
        ygl::log_fatal("at least argument should be >0");
    if (res_width < 0)
        res_width =
            (int)round(img.width() * (res_height / (float)img.height()));
    if (res_height < 0)
        res_height =
            (int)round(img.height() * (res_width / (float)img.width()));
    auto res = Image(res_width, res_height);
    ygl::resize_image(img, res);
    return res;
}

template <typename Image>
Image make_image_grid(const std::vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width() * tilex;
    auto height =
        imgs[0].height() * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = ygl::make_image(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.width() != imgs[0].width() ||
            img.height() != imgs[0].height()) {
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

ygl::image4f filter_bilateral(const ygl::image4f& img, float spatial_sigma,
    float range_sigma, const std::vector<ygl::image4f>& features,
    const std::vector<float>& features_sigma) {
    auto filtered = ygl::image4f(img.width(), img.height());
    auto width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = std::vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
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
    auto filtered = ygl::image4f(img.width(), img.height());
    auto width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
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
    // command line params
    auto parser = ygl::make_parser(argc, argv, "yimproc", "process images");
    auto command = ygl::parse_arg(parser, "command", "command to execute", ""s,
        true, {"resize", "tonemap", "bilateral"});
    auto output =
        ygl::parse_opt(parser, "--output", "-o", "output image filename", ""s);
    if (command == "resize") {
        auto width = ygl::parse_opt(
            parser, "--width", "-w", "width (-1 to maintain aspect)", -1);
        auto height = ygl::parse_opt(
            parser, "--height", "-h", "height (-1 to maintain aspect)", -1);
        auto filename =
            ygl::parse_arg(parser, "filename", "input image filename", ""s);
        // check parsing
        if (ygl::should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        if (ygl::is_hdr_filename(filename)) {
            auto img = load_hdr(filename);
            auto out = resize_image(img, width, height);
            save_hdr(output, out);
        } else {
            auto img = load_ldr(filename);
            auto out = resize_image(img, width, height);
            save_ldr(output, out);
        }
    } else if (command == "tonemap") {
        auto params = ygl::parse_params(parser, "", ygl::tonemap_params());
        auto filename =
            ygl::parse_arg(parser, "filename", "input image filename", ""s);
        // check parsing
        if (ygl::should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        auto img = load_hdr(filename);
        auto out = tonemap_image(img, params);
        save_ldr(output, out);
    } else if (command == "bilateral") {
        auto spatial_sigma =
            ygl::parse_opt(parser, "--spatial-sigma", "-s", "spatial sigma", 3);
        auto range_sigma =
            ygl::parse_opt(parser, "--range-sigma", "-r", "range sigma", 0.1f);
        auto feature_sigma = ygl::parse_opt(
            parser, "--features-sigma", "-f", "features sigmas", 0.1f);
        auto filename =
            ygl::parse_arg(parser, "filename", "input image filename", ""s);
        auto ffilenames = ygl::parse_args(parser, "features",
            "input features filename", std::vector<std::string>{}, false);
        // check parsing
        if (ygl::should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        auto img = load_hdr(filename);
        auto features = std::vector<ygl::image4f>();
        auto features_sigma = std::vector<float>();
        for (auto ffilename : ffilenames) {
            features.push_back(load_hdr(ffilename));
            features_sigma.push_back(feature_sigma);
        }
        auto out = filter_bilateral(
            img, spatial_sigma, range_sigma, features, features_sigma);
        save_hdr(output, out);
    } else {
        // check parsing
        if (ygl::should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }
    }

    // done
    return 0;
}
