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

//
// Compare two images are returns either 0 or 1.
//

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

image4f compute_diff_image(const image4f& a, const image4f& b) {
    auto diff = image{a.width, a.height, zero4f};
    for (auto j = 0; j < a.height; j++) {
        for (auto i = 0; i < a.width; i++) {
            diff[{i, j}] = {abs(a[{i, j}].x - b[{i, j}].x),
                abs(a[{i, j}].y - b[{i, j}].y), abs(a[{i, j}].z - b[{i, j}].z),
                abs(a[{i, j}].w - b[{i, j}].w)};
        }
    }
    return diff;
}

vec4f max_diff_value(const image4f& diff) {
    auto max_value = vec4f{0, 0, 0, 0};
    for (auto& c : diff) {
        max_value = {max(c.x, max_value.x), max(c.y, max_value.y),
            max(c.z, max_value.z), max(c.w, max_value.w)};
    }
    return max_value;
}

image4f display_diff(const image4f& diff) {
    auto display = image{diff.width, diff.height, zero4f};
    for (auto j = 0; j < diff.height; j++) {
        for (auto i = 0; i < diff.width; i++) {
            auto diff_value = max(diff[{i, j}]);
            display[{i, j}] = {diff_value, diff_value, diff_value, 1};
        }
    }
    return display;
}

int main(int argc, char* argv[]) {
    // parse command line
    auto parser = cmdline_parser{};
    init_cmdline_parser(parser, argc, argv, "Compares two images", "yimdiff");
    auto threshold = parse_cmdline_argument(
        parser, "--threshold,-t", 0.1f, "Thhhreshold");
    auto output = parse_cmdline_argument(
        parser, "--output,-o", ""s, "output image filename", false);
    auto filename1 = parse_cmdline_argument(
        parser, "filename1", "in1.png"s, "input image filename", true);
    auto filename2 = parse_cmdline_argument(
        parser, "filename2", "in2.png"s, "input image filename", true);
    check_cmdline_parser(parser);

    // check image type
    auto img1 = image4f{}, img2 = image4f{};
    if (!load_image(filename1, img1))
        log_fatal("cannot open image {}", filename1);
    if (!load_image(filename2, img2))
        log_fatal("cannot open image {}", filename2);
    if (img1.width != img2.width || img1.height != img2.height)
        log_fatal("image size differs");
    auto diff     = compute_diff_image(img1, img2);
    auto max_diff = max_diff_value(diff);
    if (!empty(output)) {
        auto display = display_diff(diff);
        if (!save_image(output, display))
            log_fatal("cannot save image {}", output);
    }
    if (max(max_diff) > threshold) {
        log_info("image max difference: {}", max_diff);
        log_fatal("image content differs");
    }

    // done
    return 0;
}
