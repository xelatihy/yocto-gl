//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

int main(int argc, const char* argv[]) {
  // command line parameters
  auto smooth    = false;
  auto scale     = vec3f{1, 1, 1};
  auto uscale    = 1.0f;
  auto height    = 1.0f;
  auto rotate    = zero3f;
  auto translate = zero3f;
  auto info      = false;
  auto output    = "out.ply"s;
  auto filename  = "heightfield.png"s;

  // parse command line
  auto cli = make_cli("yheightfieldproc", "Makes a mesh from a heightfield");
  add_optional(cli, "height", height, "Height scale", "h");
  add_optional(cli, "smooth", smooth, "Compute smooth normals");
  add_optional(cli, "rotatey", rotate.y, "Rotate around y axis", "ry");
  add_optional(cli, "rotatex", rotate.x, "Rotate around x axis", "rx");
  add_optional(cli, "rotatez", rotate.z, "Rotate around z axis", "rz");
  add_optional(cli, "translatey", translate.y, "Translate along y axis", "ty");
  add_optional(cli, "translatex", translate.x, "Translate along x axis", "tx");
  add_optional(cli, "translatez", translate.z, "Translate along z axis", "tz");
  add_optional(cli, "scale", uscale, "Scale along xyz axes", "s");
  add_optional(cli, "scaley", scale.y, "Scale along y axis", "sy");
  add_optional(cli, "scalex", scale.x, "Scale along x axis", "sx");
  add_optional(cli, "scalez", scale.z, "Scale along z axis", "sz");
  add_optional(cli, "info", info, "print mesh info", "i");
  add_optional(cli, "output", output, "output mesh", "o");
  add_positional(cli, "mesh", filename, "input heightfield");
  parse_cli(cli, argc, argv);

  // mesh data
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto texcoords = vector<vec2f>{};
  auto quads     = vector<vec4i>{};

  // load mesh
  auto img     = image<vec4f>{};
  auto ioerror = ""s;
  print_progress("load image", 0, 1);
  if (is_hdr_filename(filename)) {
    if (!load_image(filename, img, ioerror)) print_fatal(ioerror);
  } else {
    auto img16 = image<vec4s>{};
    if (!load_image(filename, img16, ioerror)) print_fatal(ioerror);
    img = ushort_to_float(img16);
  }
  print_progress("load image", 1, 1);

  // heightfield data
  auto heightfield = rgba_to_gray(img);

  // adjust height
  if (height != 1) {
    for (auto& pixel : heightfield) pixel *= height;
  }

  // create heightfield
  make_heightfield(quads, positions, normals, texcoords, heightfield.imsize(),
      heightfield.data_vector());
  if (!smooth) normals.clear();

  // print info
  if (info) {
    print_info("shape stats ------------");
    auto stats = shape_stats(
        {}, {}, {}, quads, {}, {}, {}, positions, normals, texcoords, {}, {});
    for (auto& stat : stats) print_info(stat);
  }

  // transform
  if (uscale != 1) scale *= uscale;
  if (translate != zero3f || rotate != zero3f || scale != vec3f{1, 1, 1}) {
    print_progress("transform shape", 0, 1);
    auto xform = translation_frame(translate) * scaling_frame(scale) *
                 rotation_frame({1, 0, 0}, radians(rotate.x)) *
                 rotation_frame({0, 0, 1}, radians(rotate.z)) *
                 rotation_frame({0, 1, 0}, radians(rotate.y));
    for (auto& p : positions) p = transform_point(xform, p);
    for (auto& n : normals)
      n = transform_normal(xform, n, max(scale) != min(scale));
    print_progress("transform shape", 1, 1);
  }

  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_shape(output, {}, {}, {}, quads, {}, {}, {}, positions, normals,
          texcoords, {}, {}, ioerror))
    print_fatal(ioerror);
  print_progress("save shape", 1, 1);

  // done
  return 0;
}
