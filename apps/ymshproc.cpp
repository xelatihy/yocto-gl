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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

void exit_error(const string& msg) {
    printf("%s\n", msg.c_str());
    exit(1);
}

int main(int argc, char** argv) {
    // command line parameters
    auto geodesic_source = -1;
    auto output          = "out.ply"s;
    auto filename        = "mesh.ply"s;

    // parse command line
    auto parser = CLI::App{"Applies operations on a triangle mesh"};
    parser.add_option(
        "--geodesic-source,-g", geodesic_source, "Geodesic source");
    parser.add_option("--output,-o", output, "output mesh")->required(true);
    parser.add_option("mesh", filename, "input mesh")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // load mesh
    auto shape = yocto_shape{};
    try {
        load_mesh(filename, shape.points, shape.lines, shape.triangles,
            shape.quads, shape.positions, shape.normals, shape.texturecoords,
            shape.colors, shape.radius, true);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // compute geodesics and store them as colors
    if (geodesic_source >= 0) {
        auto solver = geodesic_solver<float>{};
        init_geodesic_solver(solver, shape.triangles, shape.positions);
        auto distances = vector<float>{};
        compute_geodesic_distances(solver, distances, {geodesic_source});
        shape.colors = vector<vec4f>{};
        convert_distance_to_color(shape.colors, distances);
    }

    // save mesh
    try {
        save_mesh(output, shape.points, shape.lines, shape.triangles,
            shape.quads, shape.positions, shape.normals, shape.texturecoords,
            shape.colors, shape.radius);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // done
    return 0;
}
