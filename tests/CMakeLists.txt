if(YOCTO_TESTING)
  enable_testing()
  set(test-yimgproc ${CMAKE_SOURCE_DIR}/tests/_results/test-yimgproc.sh)
  set(test-yscntrace ${CMAKE_SOURCE_DIR}/tests/_results/test-yscntrace.sh)
  set(test-yscnproc ${CMAKE_SOURCE_DIR}/tests/_results/test-yscnproc.sh)
  set(testdir ${CMAKE_SOURCE_DIR}/tests)
  set(outdir ${CMAKE_SOURCE_DIR}/tests/_output)
  set(resdir ${CMAKE_SOURCE_DIR}/tests/_results)
  set(diffdir ${CMAKE_SOURCE_DIR}/tests/_difference)
  # image io
  foreach(filename "uvgrid.png")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yimgproc-${filename}
            COMMAND ${test-yimgproc} ${testdir}/textures/${filename} ${outdir}/yimgproc-${name}.png ${resdir}/yimgproc-${name}.png ${diffdir}/yimgproc-${name}.png
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})    
  endforeach()
  # procedurla images
  foreach(name "images1" "images2")
    add_test(NAME yimgproc-${name}
            COMMAND ${test-yimgproc} ${name}.ypreset ${outdir}/yimgproc-${name}.png ${resdir}/yimgproc-${name}.png ${diffdir}/yimgproc-${name}.png
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})    
  endforeach()
  # eyelight
  foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "materials4.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "eyelight")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 4 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  # path
  foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "materials4.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "path")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  # naive
  foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "naive")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  if(YOCTO_EMBREE)
  # embree
  foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "path")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-embree-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-embree-${tracer}-${name}.png ${resdir}/yscntrace-embree-${tracer}-${name}.png ${diffdir}/yscntrace-embree-${tracer}-${name}.png "-t ${tracer} --bvh-embree -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  endif(YOCTO_EMBREE)
  # yaml output
  foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "yaml")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  # obj output
  foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "obj")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
  # pbrt output
  foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "pbrt")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endforeach()
endif(YOCTO_TESTING)
