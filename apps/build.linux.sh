mkdir -p bin

cc="clang"
cflags="-Ofast -ffast-math -funroll-loops -std=gnu99"
linkflags="-lm -lpthread"

$cc -DYA_NOGL $cflags $linkflags apps/ytrace.c -o bin/ytrace-cli
$cc $cflags $linkflags apps/ytestgen.c -o bin/ytestgen

