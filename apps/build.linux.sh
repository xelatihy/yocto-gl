mkdir -p bin

g++-6 -MMD -MF bin/ysym.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ysym.o apps/ysym.cpp
g++-6 -o bin/ysym bin/ysym.o --std=c++14 -pthread
g++-6 -MMD -MF bin/ytestgen.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytestgen.o apps/ytestgen.cpp
g++-6 -o bin/ytestgen bin/ytestgen.o --std=c++14 -pthread
g++-6 -MMD -MF bin/ytrace.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytrace.o apps/ytrace.cpp
g++-6 -o bin/ytrace bin/ytrace.o --std=c++14 -pthread
g++-6 -MMD -MF bin/yobj2gltf.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yobj2gltf.o apps/yobj2gltf.cpp
g++-6 -o bin/yobj2gltf bin/yobj2gltf.o --std=c++14 -pthread
