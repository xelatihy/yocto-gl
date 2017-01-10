CXX=clang++
CXXFLAGS='-std=c++14 -stdlib=libc++ -O3 -g'
OPENGL_LFLAGS='-I/usr/local/include -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -L/usr/local/lib -lglfw3'

mkdir -p bin

$CXX $CXXFLAGS -o bin/ysym apps/ysym.cpp
$CXX $CXXFLAGS -o bin/ytestgen apps/ytestgen.cpp
$CXX $CXXFLAGS -o bin/ytrace apps/ytrace.cpp
$CXX $CXXFLAGS -o bin/yobj2gltf apps/yobj2gltf.cpp

$CXX $CXXFLAGS $OPENGL_LFLAGS -o bin/yshade apps/yshade.cpp
$CXX $CXXFLAGS $OPENGL_LFLAGS -o bin/yisym apps/yisym.cpp
$CXX $CXXFLAGS $OPENGL_LFLAGS -o bin/yitrace apps/yitrace.cpp
$CXX $CXXFLAGS $OPENGL_LFLAGS -o bin/yimview apps/yimview.cpp
