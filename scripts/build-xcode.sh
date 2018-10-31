mkdir -p build
mkdir -p build/xcode
cd build/xcode
cmake ../.. -GXcode -DYOCTO_EMBREE=ON
open yocto-gl.xcodeproj
