mkdir -p bin

g++ -MMD -MF bin/ysym.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ysym.o apps/ysym.cpp
g++ -MMD -MF bin/yapp.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yapp.o apps/yapp.cpp
g++ -MMD -MF bin/tinyply.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/tinyply.o apps/tinyply.cpp
g++ -MMD -MF bin/yocto_bvh.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_bvh.o yocto/yocto_bvh.cpp
g++ -MMD -MF bin/yocto_gltf.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_gltf.o yocto/yocto_gltf.cpp
g++ -MMD -MF bin/yocto_obj.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_obj.o yocto/yocto_obj.cpp
g++ -MMD -MF bin/yocto_trace.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_trace.o yocto/yocto_trace.cpp
g++ -MMD -MF bin/yocto_sym.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_sym.o yocto/yocto_sym.cpp
g++ -MMD -MF bin/yocto_shape.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_shape.o yocto/yocto_shape.cpp
g++ -MMD -MF bin/yocto_cmd.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yocto_cmd.o yocto/yocto_cmd.cpp
ar rcs bin/yocto.a bin/yocto_bvh.o bin/yocto_gltf.o bin/yocto_obj.o bin/yocto_trace.o bin/yocto_sym.o bin/yocto_shape.o bin/yocto_cmd.o
g++ -o bin/ysym bin/ysym.o bin/yapp.o bin/tinyply.o bin/yocto.a --std=c++14 -pthread
g++ -MMD -MF bin/ytestgen.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytestgen.o apps/ytestgen.cpp
g++ -o bin/ytestgen bin/ytestgen.o bin/yapp.o bin/tinyply.o bin/yocto.a --std=c++14 -pthread
g++ -MMD -MF bin/ytrace.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytrace.o apps/ytrace.cpp
g++ -o bin/ytrace bin/ytrace.o bin/yapp.o bin/tinyply.o bin/yocto.a --std=c++14 -pthread
g++ -MMD -MF bin/yobj2gltf.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yobj2gltf.o apps/yobj2gltf.cpp
g++ -o bin/yobj2gltf bin/yobj2gltf.o bin/yapp.o bin/tinyply.o bin/yocto.a --std=c++14 -pthread
