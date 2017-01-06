clang++ -MMD -MF bin/yshade.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -L/usr/local/lib -lglfw3 -o bin/yshade apps/yshade.cpp
clang++ -MMD -MF bin/ysym.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -o bin/ysym apps/ysym.cpp
clang++ -MMD -MF bin/yisym.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -L/usr/local/lib -lglfw3 -o bin/yisym apps/yisym.cpp
clang++ -MMD -MF bin/ytestgen.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -o bin/ytestgen apps/ytestgen.cpp
clang++ -MMD -MF bin/ytrace.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -o bin/ytrace apps/ytrace.cpp
clang++ -MMD -MF bin/yitrace.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -L/usr/local/lib -lglfw3 -o bin/yitrace apps/yitrace.cpp
clang++ -MMD -MF bin/yimview.d -std=c++14 -stdlib=libc++ -O3 -g -fcolor-diagnostics -I/usr/local/include -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -L/usr/local/lib -lglfw3 -o bin/yimview apps/yimview.cpp
