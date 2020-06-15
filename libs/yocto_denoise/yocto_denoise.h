//
// Yocto/Denoise: Denoise images using Intel Open Image Denoise.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_DENOISE_H_
#define _YOCTO_DENOISE_H_

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

namespace yocto::denoise {

using progress_callback =
    std::function<void(const std::string &message, int current, int total)>;

// Denoise a standalone image without providing any additional feature images.
// The argument 'hdr' signals wether the 'color' image is HDR or LDR
bool oidn_image_denoise(const image<vec3f> &color, bool hdr, image<vec3f> &out,
    std::string &error, progress_callback progress_cb = {});

// Denoise an image with its corresponding albedo feature image.
// The albedo feature image can be either LDR or HDR, and it helps the denoiser
// mantain texture detail
bool oidn_image_denoise(const image<vec3f> &color, bool hdr,
    const image<vec3f> &albedo, image<vec3f> &out, std::string &error,
    progress_callback progress_cb = {});

// Denoise an image with both image and normal feature images.
// The normal feature image must be an HDR image as oidn expects the normals to
// be symmetrical to zero (range [-inf,inf]).
// The normal feature image helps the denoiser mantain edge detail.
bool oidn_image_denoise(const image<vec3f> &color, bool hdr,
    const image<vec3f> &albedo, const image<vec3f> &normal, image<vec3f> &out,
    std::string &error, progress_callback progress_cb = {});

}  // namespace yocto::denoise

#endif