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
#include "CLI11.hpp"
using namespace std::literals;

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
    // command line parameters
    auto filename = "scene.json"s;
    auto output = "output.json"s;
    auto notextures = false;

    // command line params
    CLI::App parser("scene processing utility", "yscnproc");
    parser.add_flag("--notextures", notextures, "Disable textures.");
    parser.add_option("--output,-o", output, "output scene")->required(true);
    parser.add_option("scene", filename, "input scene")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // load scene
    auto scn = std::shared_ptr<ygl::scene>();
    try {
        scn = ygl::load_scene(filename, !notextures);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // make a directory if needed
    try {
        mkdir(ygl::get_dirname(output));
    } catch (const std::exception& e) {
        std::cout << "cannot create directory " << ygl::get_dirname(output)
                  << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    // save scene
    try {
        ygl::save_scene(output, scn);
    } catch (const std::exception& e) {
        std::cout << "cannot save scene " << output << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // done
    return 0;
}
