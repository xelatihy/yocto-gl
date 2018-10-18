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

//
// Compare two images are returns either 0 or 1.
//

#include "../../yocto/ygl.h"
#include "../../yocto/yglio.h"
using namespace ygl;

int main(int argc, char* argv[]) {
    // parse command line
    auto parser  = make_cmdline_parser(argc, argv, "Compares two images", "yimdiff");
    auto output = parse_arg(
        parser, "--output,-o", "out.png"s, "output image filename", false);
    auto filename1 = parse_arg(
        parser, "filename1", "in1.png"s, "input image filename", true);
    auto filename2 = parse_arg(
        parser, "filename2", "in2.png"s, "input image filename", true);
    check_cmdline(parser);

    // check image type
    if(is_hdr_filename(filename1) && is_hdr_filename(filename2)) {
        auto img1 = load_image4f(filename1);
        if(img1.pixels.empty()) log_fatal("cannot open image {}", filename1);
        auto img2 = load_image4f(filename2);
        if(img2.pixels.empty()) log_fatal("cannot open image {}", filename2);
        if(img1.width != img2.width || img1.height != img2.height) 
            log_fatal("image size differs");
        if(img1.pixels != img2.pixels) 
            log_fatal("image content differs");
        // if(!output.empty()) {
        //     if(!save_image4f(output, diff)) log_fatal("cannot save image {}", output);
        // }
    } else if(!is_hdr_filename(filename1) && !is_hdr_filename(filename2)) {
        auto img1 = load_image4b(filename1);
        if(img1.pixels.empty()) log_fatal("cannot open image {}", filename1);
        auto img2 = load_image4b(filename2);
        if(img2.pixels.empty()) log_fatal("cannot open image {}", filename2);
        if(img1.width != img2.width || img1.height != img2.height) 
            log_fatal("image size differs");
        if(img1.pixels != img2.pixels) 
            log_fatal("image content differs");
        // if(!output.empty()) {
        //     if(!save_image4f(output, diff)) log_fatal("cannot save image {}", output);
        // }
    } else {
        log_fatal("different image types");
    }

    // done
    return 0;
}
