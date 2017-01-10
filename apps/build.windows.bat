mkdir bin

set CXX=cl
set CXXFLAGS=/Od /Zi /Fdbin\ /Fobin\ /Febin\ /EHsc /bigobj
set OPENGL_LFLAGS=/Iapps/w32/include opengl32.lib apps/w32/lib-vc2015/glew32.lib apps/w32/lib-vc2015/glfw3dll.lib

%CXX% %CXXFLAGS% apps/ysym.cpp
%CXX% %CXXFLAGS% apps/ytrace.cpp
%CXX% %CXXFLAGS% apps/ytestgen.cpp
%CXX% %CXXFLAGS% apps/yobj2gltf.cpp

%CXX% %CXXFLAGS% %OPENGL_LFLAGS% apps/yshade.cpp
%CXX% %CXXFLAGS% %OPENGL_LFLAGS% apps/yisym.cpp
%CXX% %CXXFLAGS% %OPENGL_LFLAGS% apps/yitrace.cpp
%CXX% %CXXFLAGS% %OPENGL_LFLAGS% apps/yimview.cpp
