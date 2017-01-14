mkdir -p bin

clang++ -MMD -MF bin/ysym.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ysym.o apps/ysym.cpp
clang++ -MMD -MF bin/yapp.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yapp.o apps/yapp.cpp
clang++ -MMD -MF bin/yocto_bvh.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_bvh.o yocto/yocto_bvh.cpp
clang++ -MMD -MF bin/yocto_gltf.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_gltf.o yocto/yocto_gltf.cpp
clang++ -MMD -MF bin/yocto_obj.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_obj.o yocto/yocto_obj.cpp
clang++ -MMD -MF bin/yocto_trace.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_trace.o yocto/yocto_trace.cpp
clang++ -MMD -MF bin/yocto_sym.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_sym.o yocto/yocto_sym.cpp
clang++ -MMD -MF bin/yocto_shape.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yocto_shape.o yocto/yocto_shape.cpp
libtool -static bin/yocto_bvh.o bin/yocto_gltf.o bin/yocto_obj.o bin/yocto_trace.o bin/yocto_sym.o bin/yocto_shape.o -o bin/yocto.a
clang++ -o bin/ysym bin/ysym.o bin/yapp.o bin/yocto.a -std=c++14
clang++ -MMD -MF bin/ytestgen.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ytestgen.o apps/ytestgen.cpp
clang++ -o bin/ytestgen bin/ytestgen.o bin/yapp.o bin/yocto.a -std=c++14
clang++ -MMD -MF bin/ytrace.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/ytrace.o apps/ytrace.cpp
clang++ -o bin/ytrace bin/ytrace.o bin/yapp.o bin/yocto.a -std=c++14
clang++ -MMD -MF bin/yobj2gltf.o.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -c -o bin/yobj2gltf.o apps/yobj2gltf.cpp
clang++ -o bin/yobj2gltf bin/yobj2gltf.o bin/yapp.o bin/yocto.a -std=c++14
