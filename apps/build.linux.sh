mkdir -p bin

g++ -MMD -MF bin/ysym.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ysym.o apps/ysym.cpp
g++ -o bin/ysym bin/ysym.o --std=c++14 -pthread
g++ -MMD -MF bin/ytestgen.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytestgen.o apps/ytestgen.cpp
g++ -o bin/ytestgen bin/ytestgen.o --std=c++14 -pthread
g++ -MMD -MF bin/ytrace.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/ytrace.o apps/ytrace.cpp
g++ -o bin/ytrace bin/ytrace.o --std=c++14 -pthread
g++ -MMD -MF bin/yobj2gltf.o.d -std=c++14 -O3 -g -pthread -I/usr/include -c -o bin/yobj2gltf.o apps/yobj2gltf.cpp
g++ -o bin/yobj2gltf bin/yobj2gltf.o --std=c++14 -pthread
