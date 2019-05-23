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
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

image<vec4f> compute_diff_image(const image<vec4f>& a, const image<vec4f>& b) {
    auto diff = image{a.size(), zero4f};
    for (auto j = 0; j < a.size().y; j++) {
        for (auto i = 0; i < a.size().x; i++) {
            diff[{i, j}] = {abs(a[{i, j}].x - b[{i, j}].x),
                abs(a[{i, j}].y - b[{i, j}].y), abs(a[{i, j}].z - b[{i, j}].z),
                abs(a[{i, j}].w - b[{i, j}].w)};
        }
    }
    return diff;
}

vec4f max_diff_value(const image<vec4f>& diff) {
    auto max_value = vec4f{0, 0, 0, 0};
    for (auto& c : diff) {
        max_value = {max(c.x, max_value.x), max(c.y, max_value.y),
            max(c.z, max_value.z), max(c.w, max_value.w)};
    }
    return max_value;
}

image<vec4f> display_diff(const image<vec4f>& diff) {
    auto display = image{diff.size(), zero4f};
    for (auto j = 0; j < diff.size().y; j++) {
        for (auto i = 0; i < diff.size().x; i++) {
            auto diff_value = max(diff[{i, j}]);
            display[{i, j}] = {diff_value, diff_value, diff_value, 1};
        }
    }
    return display;
}

int main(int argc, char* argv[]) {
    // Application values
    auto threshold = 0.1f;
    auto output    = ""s;
    auto filename1 = ""s;
    auto filename2 = ""s;

    // parse command line
    auto parser = CLI::App{"Compare two images"};
    parser.add_option("--threshold,-t", threshold, "Threshold");
    parser.add_option("--output,-o", output, "output image filename");
    parser.add_option("filename1", filename1, "input image filename")
        ->required(true);
    parser.add_option("filename2", filename2, "input image filename")
        ->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // check image type
    auto img1 = image<vec4f>{}, img2 = image<vec4f>{};
    try {
        load_image(filename1, img1);
        load_image(filename2, img2);
    } catch (const std::exception& e) {
        print_fatal(e.what());
    }
    if (img1.size() != img2.size()) print_fatal("image size differs");
    auto diff     = compute_diff_image(img1, img2);
    auto max_diff = max_diff_value(diff);
    if (!output.empty()) {
        auto display = display_diff(diff);
        try {
            save_image(output, display);
        } catch (const std::exception& e) {
            print_fatal(e.what());
        }
    }
    if (max(max_diff) > threshold) {
        print_info("image max difference: " + to_string(max_diff));
        print_fatal("image content differs");
    }

    // done
    return 0;
}
