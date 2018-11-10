mkdir -p build
mkdir -p build/release
cd build/release
cmake ../.. -GNinja -DYOCTO_EMBREE=ON
cmake --build .
rm tests/run-tests/output/*.png
rm tests/run-tests/difference/*.png
ctest -j 4 --output-on-failure $1
