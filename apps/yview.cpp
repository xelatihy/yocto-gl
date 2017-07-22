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

// ---------------------------------------------------------------------------
// SCENE
// ---------------------------------------------------------------------------

#include "yscene.h"

// ---------------------------------------------------------------------------
// INTERACTIVE FUNCTIONS
// ---------------------------------------------------------------------------

bool update(yscene* scn) {
    // advance time
    if (scn->animate && scn->time_range != ym::zero2f) {
        if (scn->time >= scn->time_range.y) {
            scn->time = scn->time_range.x;
        } else {
            scn->time += 1 / 60.0f;
        }
        return true;
    }
    return false;
}

void draw_custom_widgets(ygui::window* win) {
    auto scn = (yscene*)ygui::get_user_pointer(win);
    if (scn->time_range != ym::zero2f) {
        ygui::slider_widget(
            win, "time", &scn->time, scn->time_range.x, scn->time_range.y);
        ygui::checkbox_widget(win, "play animation", &scn->animate);
    }
}

// ---------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // create empty scene
    auto scn = new yscene();

    // params
    parse_cmdline(
        scn, argc, argv, "yview", "interactively view scenes", false, false);

    // scene
    if (!load_scene(scn, scn->filename, false, false)) return 1;

    // run ui
    auto width = (int)std::round(scn->view_cam->aspect * scn->resolution);
    auto height = scn->resolution;
    run_ui(scn, width, height, "yview", shade_init, shade_draw, update);

    // clear
    delete scn;

    // done
    return 0;
}
