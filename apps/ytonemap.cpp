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

int main(int argc, char* argv[]) {
    // parse command line
    auto parser =
        ygl::make_cmdline_parser(argc, argv, "Process images", "ytonemap");
    auto exposure =
        ygl::parse_float(parser, "--exposure,-e", 0, "Tonemap exposure");
    auto gamma = ygl::parse_float(parser, "--gamma,-g", 2.2f, "Tonemap gamma.");
    auto filmic = ygl::parse_flag(
        parser, "--filmic,-f", false, "Tonemap uses filmic curve");
    auto output = ygl::parse_string(
        parser, "--output,-o", "out.png", "output image filename", true);
    auto filename = ygl::parse_string(
        parser, "filename", "img.hdr", "input image filename", true);
    ygl::check_cmdline(parser);

    auto hdr = ygl::load_image4f(filename);
    auto ldr = ygl::tonemap_image4f(hdr, exposure, gamma, filmic);
    ygl::save_image4f(output, ldr);

    // done
    return 0;
}
