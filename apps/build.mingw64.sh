mkdir -p bin

cc="gcc"
cflags="-Ofast -ffast-math -funroll-loops -std=gnu99"
linkflags="-lglfw3 -lOpenGL32 -lglew32 -lpthread"

# at present, ytrace fails to compile on windows without MSVC
#$cc $cflags apps/ytrace.c $linkflags -o bin/ytrace-cli

$cc $cflags apps/ytestgen.c $linkflags -o bin/ytestgen
$cc $cflags apps/yshade.c $linkflags -o bin/yshade
