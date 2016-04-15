mkdir -p bin

brew_prefix=`brew --prefix`
brew_prefix=${brew_prefix:-"/Users/fabio/homebrew"}

cc="clang"
cflags="-I$brew_prefix/include/ -Ofast -ffast-math -funroll-loops -march=native -fdiagnostics-color=always -Wall -Wpedantic"
linkflags="-L$brew_prefix/lib -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo"

cp apps/render_tests.sh bin/

$cc $cflags $linkflags apps/yview.c -o bin/yview &
$cc $cflags $linkflags apps/ytrace.c -o bin/ytrace &
$cc -DYA_NOGL $cflags $linkflags apps/ytrace.c -o bin/ytrace-cli &
$cc $cflags $linkflags apps/yshade.c -o bin/yshade &
$cc $cflags $linkflags apps/ytestgen.c -o bin/ytestgen &
wait
