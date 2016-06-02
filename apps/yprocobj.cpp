//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"

yo_scene* load(const char* filename) {
    char ext[16];
    yc_split_path(filename, 0, 0, ext);
    if (!strcmp(ext, ".obj")) return yo_load_obj(filename, false, true);
    if (!strcmp(ext, ".objbin")) return yo_load_objbin(filename, true);
    assert(false);
    return nullptr;
}

bool save(const char* filename, const yo_scene* scene) {
    char ext[16];
    yc_split_path(filename, 0, 0, ext);
    if (!strcmp(ext, ".obj")) return yo_save_obj(filename, scene, true);
    if (!strcmp(ext, ".objbin")) return yo_save_objbin(filename, scene, true);
    assert(false);
    return false;
}

int main(int argc, const char** argv) {
    // command line params
    yc_parser* parser = yc_init_parser(argc, argv, "make tests");
    const char* filename_in =
        yc_parse_args(parser, "filename_in", "input filename", 0, true);
    const char* filename_out =
        yc_parse_args(parser, "filename_out", "output filename", 0, true);
    yc_done_parser(parser);

    yo_scene* scene = load(filename_in);

    save(filename_out, scene);
}
