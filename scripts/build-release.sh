mkdir -p build
mkdir -p build/release
cd build/release
cmake ../.. -GNinja -DYOCTO_EMBREE=ON
cmake --build .
