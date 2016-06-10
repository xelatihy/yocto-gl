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

// #define YGL_USESTL

#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_obj.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "ext/tiny_obj_loader.h"

int main(int argc, const char** argv) {
    // command line
    yc_parser* parser = yc_init_parser(argc, argv, "test obj");
    const char* filename =
        yc_parse_args(parser, "scene", "scene filename", 0, true);
    yc_done_parser(parser);

    // load scene
    ym_timer yocto_timer = ym_timer();
    yo_scene* scene = yo_load_obj(filename, true, true);
    double yocto_elapsed = yocto_timer.elapsed();
    printf("yocto obj: %d shapes in %lg ms\n", (int)scene->shapes.size(),
           yocto_elapsed);
    yo_free_scene(scene);

    // tinyobj load scene
    ym_timer tiny_timer = ym_timer();
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    tinyobj::LoadObj(shapes, materials, err, filename);
    double tiny_elapsed = tiny_timer.elapsed();
    printf("tiny  obj: %d shapes in %lg ms\n", (int)shapes.size(),
           tiny_elapsed);

    // done
    return EXIT_SUCCESS;
}
