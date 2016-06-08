#define YGL_DECLARATION
#define YGL_IMPLEMENTATION

// #define YGL_USESTL

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include "ext/glew/glew.h"
#endif

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_symrigid.h"
#include "../yocto/yocto_trace.h"
