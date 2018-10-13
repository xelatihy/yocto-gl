mkdir -p build
mkdir -p build/xcode
cd build/xcode
cmake ../.. -GXcode
open yocto-gl.xcodeproj
