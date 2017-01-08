mkdir -p bin
# g++ -MMD -MF bin/yshade.d -std=c++14 -O3 -g -lgl -lglew -lglfw3 -o bin/yshade apps/yshade.cpp
g++ -MMD -MF bin/ysym.d -std=c++14 -O3 -g -o bin/ysym apps/ysym.cpp
# g++ -MMD -MF bin/yisym.d -std=c++14 -O3 -g -lgl -lglew -lglfw3 -o bin/yisym apps/yisym.cpp
g++ -MMD -MF bin/ytestgen.d -std=c++14 -O3 -g -o bin/ytestgen apps/ytestgen.cpp
g++ -MMD -MF bin/ytrace.d -std=c++14 -O3 -g -pthread -o bin/ytrace apps/ytrace.cpp
# g++ -MMD -MF bin/yitrace.d -std=c++14 -O3 -g -lgl -lglew -lglfw3 -o bin/yitrace apps/yitrace.cpp
# g++ -MMD -MF bin/yimview.d -std=c++14 -O3 -g -lgl -lglew -lglfw3 -o bin/yimview apps/yimview.cpp
g++ -MMD -MF bin/yobj2gltf.d -std=c++14 -O3 -g -o bin/yobj2gltf apps/yobj2gltf.cpp
