mkdir -p bin

clang++ -MMD -MF bin/ysym.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ysym.o apps/ysym.cpp
clang++ -o bin/ysym bin/ysym.o -std=c++14
clang++ -MMD -MF bin/ytestgen.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ytestgen.o apps/ytestgen.cpp
clang++ -o bin/ytestgen bin/ytestgen.o -std=c++14
clang++ -MMD -MF bin/ytrace.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ytrace.o apps/ytrace.cpp
clang++ -o bin/ytrace bin/ytrace.o -std=c++14
clang++ -MMD -MF bin/yobj2gltf.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yobj2gltf.o apps/yobj2gltf.cpp
clang++ -o bin/yobj2gltf bin/yobj2gltf.o -std=c++14
