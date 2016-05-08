#!/bin/bash
mkdir -p bin

cc="clang"
cflags="-Ofast -ffast-math -funroll-loops -std=gnu99"
linkflags="-lm -lpthread -lglfw3 -lGLEW -lGL -lX11 -lXrandr -lXinerama -lXi -lXxf86vm -lXcursor -ldl"

$cc $cflags apps/yview.c $linkflags -o bin/yview
$cc $cflags apps/ytrace.c $linkflags -o bin/ytrace
$cc -DYA_NOGL $cflags apps/ytrace.c $linkflags -o bin/ytrace-cli
$cc $cflags apps/yshade.c $linkflags -o bin/yshade
$cc $cflags apps/ysym.c $linkflags -o bin/ysym
$cc $cflags apps/ytestgen.c $linkflags -o bin/ytestgen
