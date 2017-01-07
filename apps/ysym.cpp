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

#include "ysym.h"

int main(int argc, char* argv[]) {
    // command line
    auto pars = new ysym_app::params();
    auto parser = ycmd::make_parser(argc, argv, "rigid body simulation");
    ysym_app::init_params(pars, parser);
    ycmd::check_parser(parser);

    // simulate each frame and save the results to a new scene
    printf("rigid body simulation for %s to %s\n", pars->filename.c_str(),
           pars->outfilename.c_str());
    printf("simulating ...");
    fflush(stdout);
    for (auto i = 0; i < pars->nframes; i++) {
        printf("\rsimulating frame %d/%d", i, pars->nframes);
        fflush(stdout);
        ysym_app::simulate_step(pars->scene, pars->rigid_scene, pars->dt);
        std::string errmsg;
        char frame_filename[4096];
        sprintf(frame_filename, pars->outfilename.c_str(), i);
        yapp::save_scene(frame_filename, pars->scene);
    }
    printf("\rsimulating done\n");
    fflush(stdout);

    // done
    delete pars;
    return 0;
}
