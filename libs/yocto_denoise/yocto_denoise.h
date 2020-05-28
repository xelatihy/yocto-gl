#ifndef _YOCTO_DENOISE_H_
#define _YOCTO_DENOISE_H_

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

namespace yocto::denoise {

// namespace aliases
namespace img  = yocto::image;
namespace math = yocto::math;

using img::image;
using math::vec3f;

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