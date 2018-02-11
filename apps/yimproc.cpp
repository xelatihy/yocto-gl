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
using namespace ygl;

// Load hdr
image4f load_hdr(const string& filename) {
    auto img = load_image4f(filename);
    if (!img) log_fatal("cannot load image {}", filename);
    return img;
}

// Load ldr
image4b load_ldr(const string& filename) {
    auto img = load_image4b(filename);
    if (!img) log_fatal("cannot load image {}", filename);
    return img;
}

// Save hdr
void save_hdr(const string& filename, const image4f& img) {
    if (!save_image4f(filename, img))
        log_fatal("cannot save image {}", filename);
}

// Save ldr
void save_ldr(const string& filename, const image4b& img) {
    if (!save_image4b(filename, img))
        log_fatal("cannot save image {}", filename);
}

// Resize image.
template <typename Image>
Image resize_image(const Image& img, int res_width, int res_height) {
    if (res_width < 0 && res_height < 0)
        log_fatal("at least argument should be >0");
    if (res_width < 0)
        res_width =
            (int)round(img.width() * (res_height / (float)img.height()));
    if (res_height < 0)
        res_height =
            (int)round(img.height() * (res_width / (float)img.width()));
    auto res = Image(res_width, res_height);
    resize_image(img, res);
    return res;
}

template <typename Image>
Image make_image_grid(const vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width() * tilex;
    auto height =
        imgs[0].height() * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = make_image(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.width() != imgs[0].width() ||
            img.height() != imgs[0].height()) {
            log_fatal("images of different sizes are not accepted");
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

image4f filter_bilateral(const image4f& img, float spatial_sigma,
    float range_sigma, const vector<image4f>& features,
    const vector<float>& features_sigma) {
    auto filtered = image4f(img.width(), img.height());
    auto width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
                    auto uv = vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img.at(i, j) - img.at(ii, jj);
                    auto w = exp(-dot(uv, uv) * sw) * exp(-dot(rgb, rgb) * rw);
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

image4f filter_bilateral(
    const image4f& img, float spatial_sigma, float range_sigma) {
    auto filtered = image4f(img.width(), img.height());
    auto width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto av = zero4f;
            auto aw = 0.0f;
            for (auto fj = -width; fj <= width; fj++) {
                for (auto fi = -width; fi <= width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= img.width() || jj >= img.height()) continue;
                    auto uv = vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img.at(i, j) - img.at(ii, jj);
                    auto w = exp(-dot(uv, uv) * sw) * exp(-dot(rgb, rgb) * rw);
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
    auto parser = make_parser(argc, argv, "yimproc", "process images");
    auto command = parse_arg(parser, "command", "command to execute", ""s, true,
        {"resize", "tonemap", "bilateral"});
    auto output =
        parse_opt(parser, "--output", "-o", "output image filename", ""s);
    if (command == "resize") {
        auto width = parse_opt(
            parser, "--width", "-w", "width (-1 to maintain aspect)", -1);
        auto height = parse_opt(
            parser, "--height", "-h", "height (-1 to maintain aspect)", -1);
        auto filename =
            parse_arg(parser, "filename", "input image filename", ""s);
        // check parsing
        if (should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        if (is_hdr_filename(filename)) {
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
            parse_opt(parser, "--exposure", "-e", "hdr exposure", 0.0f);
        auto gamma = parse_opt(parser, "--gamma", "-g", "hdr gamma", 2.2f);
        auto filmic = parse_flag(parser, "--filmic", "-F", "hdr filmic");
        auto filename =
            parse_arg(parser, "filename", "input image filename", ""s);
        // check parsing
        if (should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        auto img = load_hdr(filename);
        auto out = tonemap_image(img, exposure, gamma, filmic);
        save_ldr(output, out);
    } else if (command == "bilateral") {
        auto spatial_sigma =
            parse_opt(parser, "--spatial-sigma", "-s", "spatial sigma", 3);
        auto range_sigma =
            parse_opt(parser, "--range-sigma", "-r", "range sigma", 0.1f);
        auto feature_sigma = parse_opt(
            parser, "--features-sigma", "-f", "features sigmas", 0.1f);
        auto filename =
            parse_arg(parser, "filename", "input image filename", ""s);
        auto ffilenames = parse_args(parser, "features",
            "input features filename", vector<string>{}, false);
        // check parsing
        if (should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }

        auto img = load_hdr(filename);
        auto features = vector<image4f>();
        auto features_sigma = vector<float>();
        for (auto ffilename : ffilenames) {
            features.push_back(load_hdr(ffilename));
            features_sigma.push_back(feature_sigma);
        }
        auto out = filter_bilateral(
            img, spatial_sigma, range_sigma, features, features_sigma);
        save_hdr(output, out);
    } else {
        // check parsing
        if (should_exit(parser)) {
            printf("%s\n", get_usage(parser).c_str());
            exit(1);
        }
    }

    // done
    return 0;
}
