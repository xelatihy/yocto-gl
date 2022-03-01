## Copyright 2021 Jefferson Amstutz
## SPDX-License-Identifier: Apache-2.0

# NOTE: This script is only to be invoked by the EmbedPTX() function.

unset(file_contents)

foreach(obj ${OBJECTS})
  get_filename_component(obj_ext ${obj} EXT)
  get_filename_component(obj_dir ${obj} DIRECTORY)

  list(POP_FRONT SYMBOL_NAMES obj_name)
  if(obj_ext MATCHES ".ptx")
    set(args --name ${obj_name} ${obj})
    execute_process(
      COMMAND "${BIN_TO_C_COMMAND}" ${args}
      WORKING_DIRECTORY ${obj_dir}
      RESULT_VARIABLE result
      OUTPUT_VARIABLE output
      ERROR_VARIABLE error_var
    )
    set(file_contents "${file_contents} \n${output}")
  endif()
endforeach()

file(WRITE "${OUTPUT}" "${file_contents}")
