mkdir bin

cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ysym.obj.pdb /Fobin/ysym.obj apps/ysym.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_bvh.obj.pdb /Fobin/yocto_bvh.obj yocto/yocto_bvh.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_gltf.obj.pdb /Fobin/yocto_gltf.obj yocto/yocto_gltf.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_obj.obj.pdb /Fobin/yocto_obj.obj yocto/yocto_obj.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_trace.obj.pdb /Fobin/yocto_trace.obj yocto/yocto_trace.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_sym.obj.pdb /Fobin/yocto_sym.obj yocto/yocto_sym.cpp
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yocto_shape.obj.pdb /Fobin/yocto_shape.obj yocto/yocto_shape.cpp
lib bin/yocto_bvh.obj bin/yocto_gltf.obj bin/yocto_obj.obj bin/yocto_trace.obj bin/yocto_sym.obj bin/yocto_shape.obj /OUT:bin/yocto.lib
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yapp.obj.pdb /Fobin/yapp.obj apps/yapp.cpp
cl /Febin/ysym.exe /Fdbin/ysym.exe.pdb bin/ysym.obj bin/yocto.lib bin/yapp.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ytestgen.obj.pdb /Fobin/ytestgen.obj apps/ytestgen.cpp
cl /Febin/ytestgen.exe /Fdbin/ytestgen.exe.pdb bin/ytestgen.obj bin/yocto.lib bin/yapp.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/ytrace.obj.pdb /Fobin/ytrace.obj apps/ytrace.cpp
cl /Febin/ytrace.exe /Fdbin/ytrace.exe.pdb bin/ytrace.obj bin/yocto.lib bin/yapp.obj /Zi
cl  /O2 /Zi /EHsc /bigobj /Iapps/w32/include /c /Fdbin/yobj2gltf.obj.pdb /Fobin/yobj2gltf.obj apps/yobj2gltf.cpp
cl /Febin/yobj2gltf.exe /Fdbin/yobj2gltf.exe.pdb bin/yobj2gltf.obj bin/yocto.lib /Zi
