mkdir bin

cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ysym.obj.pdb /Fobin/ysym.obj apps/ysym.cpp
cl /Febin/ysym.exe /Fdbin/ysym.exe.pdb bin/ysym.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ytestgen.obj.pdb /Fobin/ytestgen.obj apps/ytestgen.cpp
cl /Febin/ytestgen.exe /Fdbin/ytestgen.exe.pdb bin/ytestgen.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ytrace.obj.pdb /Fobin/ytrace.obj apps/ytrace.cpp
cl /Febin/ytrace.exe /Fdbin/ytrace.exe.pdb bin/ytrace.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yobj2gltf.obj.pdb /Fobin/yobj2gltf.obj apps/yobj2gltf.cpp
cl /Febin/yobj2gltf.exe /Fdbin/yobj2gltf.exe.pdb bin/yobj2gltf.obj /Zi
