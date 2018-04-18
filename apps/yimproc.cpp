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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_utils.h"
using namespace std::literals;

// Resize image.
template <typename T>
std::vector<T> resize_image(int width, int height, const std::vector<T>& img, 
    int res_width, int res_height) {
    if (!res_width && !res_height)
        ygl::log_fatal("at least argument should be >0");
    if (!res_width)
        res_width = (int)round(width * (res_height / (float)height));
    if (!res_height)
        res_height = (int)round(height * (res_width / (float)width));
    auto res = std::vector<T>(res_width * res_height);
    ygl::resize_image(width, height, img, res_width, res_height, res);
    return res;
}

#if 0
template <typename Image>
Image make_image_grid(const std::vector<Image>& imgs, int tilex) {
    auto nimgs = (int)imgs.size();
    auto width = imgs[0].width * tilex;
    auto height = imgs[0].height * (nimgs / tilex + ((nimgs % tilex) ? 1 : 0));
    auto ret = ygl::make_image4f(width, height, (bool)imgs[0].hdr);
    auto img_idx = 0;
    for (auto& img : imgs) {
        if (img.width != imgs[0].width || img.height != imgs[0].height) {
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

std::vector<ygl::vec4f> filter_bilateral(int width, int height,
    const std::vector<ygl::vec4f>& img, float spatial_sigma, float range_sigma,
    const std::vector<std::vector<ygl::vec4f>>& features,
    const std::vector<float>& features_sigma) {
    auto filtered = std::vector<ygl::vec4f>(width * height);
    auto filter_width = (int)ceil(2.57f * spatial_sigma);
    auto sw = 1 / (2.0f * spatial_sigma * spatial_sigma);
    auto rw = 1 / (2.0f * range_sigma * range_sigma);
    auto fw = std::vector<float>();
    for (auto feature_sigma : features_sigma)
        fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto av = ygl::zero4f;
            auto aw = 0.0f;
            for (auto fj = -filter_width; fj <= filter_width; fj++) {
                for (auto fi = -filter_width; fi <= filter_width; fi++) {
                    auto ii = i + fi, jj = j + fj;
                    if (ii < 0 || jj < 0) continue;
                    if (ii >= width || jj >= height) continue;
                    auto uv = ygl::vec2f{float(i - ii), float(j - jj)};
                    auto rgb = img[i + j * width] - img[ii + width * jj];
                    auto w = (float)std::exp(-dot(uv, uv) * sw) *
                             (float)std::exp(-dot(rgb, rgb) * rw);
                    for (auto fi = 0; fi < features.size(); fi++) {
                        auto feat =
                            features[fi][i + width * j] - features[fi][ii + width * jj];
                        w *= exp(-dot(feat, feat) * fw[fi]);
                    }
                    av += w * img[ii + width * jj];
                    aw += w;
                }
            }
            filtered[i + j * width] = av / aw;
        }
    }
    return filtered;
}

std::vector<ygl::vec4f> filter_bilateral(int width, int height,
    const std::vector<ygl::vec4f>& img, float spatial_sigma, float range_sigma) {
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

void load_image(const std::string& filename, int& width, int& height,
    std::vector<ygl::vec4b>& ldr, std::vector<ygl::vec4f>& hdr,
    ygl::tonemap_type tonemapper, float exposure) {
    if (ygl::is_hdr_filename(filename)) {
        hdr = ygl::load_image4f(filename, width, height);
        ldr = ygl::tonemap_image(hdr, tonemapper, exposure);
    } else {
        ldr = ygl::load_image4b(filename, width, height);
        hdr = srgb_to_linear(ldr);
    }
    if (ldr.empty() || hdr.empty())
        ygl::log_fatal("cannot load image {}", filename);
}

void save_image(const std::string& filename, int width, int height,
    const std::vector<ygl::vec4b>& ldr, const std::vector<ygl::vec4f>& hdr) {
    if (ygl::is_hdr_filename(filename)) {
        if (!save_image4f(filename, width, height, hdr))
            ygl::log_fatal("cannot save image {}", filename);
    } else {
        if (!save_image4b(filename, width, height, ldr))
            ygl::log_fatal("cannot save image {}", filename);
    }
}

int main(int argc, char* argv[]) {
    static auto tonemap_names = std::map<ygl::tonemap_type, std::string>{
        {ygl::tonemap_type::linear, "linear"},
        {ygl::tonemap_type::gamma, "gamma"},
        {ygl::tonemap_type::srgb, "srgb"},
        {ygl::tonemap_type::filmic1, "filmic1"},
        {ygl::tonemap_type::filmic2, "filmic2"},
        {ygl::tonemap_type::filmic3, "filmic3"},
    };

    // command line params
    auto parser = ygl::make_parser(argc, argv, "yimproc", "process images");
    auto tonemapper = ygl::parse_opt(parser, "--tonemapper", "t",
        "Tonemapper type.", tonemap_names, ygl::tonemap_type::gamma);
    auto exposure =
        ygl::parse_opt(parser, "--exposure", "t", "Hdr exposure", 0.0f);
    auto resize_width = ygl::parse_opt(
        parser, "--resize-width", "-w", "width (0 to maintain aspect)", 0);
    auto resize_height = ygl::parse_opt(
        parser, "--resize-height", "-h", "height (0 to maintain aspect)", 0);
    auto multiply_color = ygl::parse_opt(parser, "--multiply-color", "",
        "multiply b y this color", ygl::vec4f{1, 1, 1, 1});
    auto spatial_sigma = ygl::parse_opt(
        parser, "--spatial-sigma", "-s", "blur spatial sigma", 0.0f);
    auto range_sigma = ygl::parse_opt(
        parser, "--range-sigma", "-r", "bilateral blur range sigma", 0.0f);
    auto set_alpha_filename = ygl::parse_opt(
        parser, "--set-alpha", "", "set alpha as this image alpha", ""s);
    auto set_color_as_alpha_filename = ygl::parse_opt(parser,
        "--set-color-as-alpha", "", "set alpha as this image color", ""s);
    auto output = ygl::parse_opt(
        parser, "--output", "-o", "output image filename", ""s, true);
    auto filename =
        ygl::parse_arg(parser, "filename", "input image filename", ""s);
    // check parsing
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load
    auto width = 0, height = 0;
    auto hdr = std::vector<ygl::vec4f>();
    auto ldr = std::vector<ygl::vec4b>();
    load_image(filename, width, height, ldr, hdr, tonemapper, exposure);

    // set alpha
    if (set_alpha_filename != "") {
        auto alpha_hdr = std::vector<ygl::vec4f>();
        auto alpha_ldr = std::vector<ygl::vec4b>();
        auto alpha_width = 0, alpha_height = 0;
        load_image(set_alpha_filename, alpha_width, alpha_height, alpha_ldr,
            alpha_hdr, tonemapper, exposure);
        if (width != alpha_width || height != alpha_height)
            ygl::log_fatal("bad image size");
        for (auto i = 0; i < ldr.size(); i++) ldr[i].w = alpha_ldr[i].w;
        for (auto i = 0; i < hdr.size(); i++) ldr[i].w = alpha_ldr[i].w;
    }

    // set alpha
    if (set_color_as_alpha_filename != "") {
        auto alpha_hdr = std::vector<ygl::vec4f>();
        auto alpha_ldr = std::vector<ygl::vec4b>();
        auto alpha_width = 0, alpha_height = 0;
        load_image(set_color_as_alpha_filename, alpha_width, alpha_height,
            alpha_ldr, alpha_hdr, tonemapper, exposure);
        if (width != alpha_width || height != alpha_height)
            ygl::log_fatal("bad image size");
        for (auto i = 0; i < ldr.size(); i++) {
            auto& p = alpha_ldr[i];
            ldr[i].w = (uint8_t)ygl::clamp(
                ((int)p.x + (int)p.y + (int)p.z) / 3, 0, 255);
        }
        for (auto i = 0; i < hdr.size(); i++) {
            auto& p = alpha_hdr[i];
            hdr[i].w = (p.x + p.y + p.z) / 3;
        }
    }

    // multiply
    if (multiply_color != ygl::vec4f{1, 1, 1, 1}) {
        for (auto& c : hdr) c *= multiply_color;
        for (auto& c : ldr)
            c = linear_to_srgb(srgb_to_linear(c) * multiply_color);
    }

    // resize
    if (resize_width || resize_height) {
        hdr = resize_image(width, height, hdr, resize_width, resize_height);
        ldr = resize_image(width, height, ldr, resize_width, resize_height);
    }

    // bilateral
    if (spatial_sigma && range_sigma) {
        hdr = filter_bilateral(
            width, height, hdr, spatial_sigma, range_sigma, {}, {});
        ldr = tonemap_image(hdr, tonemapper, exposure);
    }

    // save
    save_image(output, width, height, ldr, hdr);

    // done
    return 0;
}
