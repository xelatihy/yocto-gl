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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "CLI11.hpp"
using namespace std::literals;

struct app_state {
    std::string filename = "scene.json"s;
    std::string output = "output.json"s;
    bool notextures = false;

    std::shared_ptr<ygl::scene> scn = nullptr;
};

void mkdir(const std::string& dir) {
    if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
        return;
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
#else
    system(("mkdir " + dir).c_str());
#endif
}

int main(int argc, char** argv) {
    auto app = std::make_shared<app_state>();

    // command line params
    CLI::App parser("scene processing utility", "yscnproc");
    parser.add_flag("--notextures", app->notextures, "Disable textures.");
    parser.add_option("--output,-o", app->output, "output scene")
        ->required(true);
    parser.add_option("scene", app->filename, "input scene")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // load obj
    try {
        app->scn = ygl::load_scene(app->filename, !app->notextures);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << app->filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // make a directory if needed
    try {
        mkdir(ygl::get_dirname(app->output));
    } catch (const std::exception& e) {
        std::cout << "cannot create directory " << ygl::get_dirname(app->output)
                  << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    // save scene
    try {
        ygl::save_scene(app->output, app->scn);
    } catch (const std::exception& e) {
        std::cout << "cannot save scene " << app->output << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // done
    return 0;
}
