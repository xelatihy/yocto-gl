mkdir bin

cl /DYA_NOGL /Ox /Fobin/ytrace.obj /Febin/ytrace-cli.exe apps/ytrace.c
cl /Ox /Fobin/ytestgen.obj /Febin/ytestgen apps/ytestgen.c
