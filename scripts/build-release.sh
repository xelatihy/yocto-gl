mkdir -p build
mkdir -p build/release
cd build/release
cmake ../.. -GNinja
cmake --build .
