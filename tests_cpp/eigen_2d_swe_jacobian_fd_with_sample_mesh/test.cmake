include(FindUnixCommands)

#------------------
# make full mesh
#------------------
set(CMD "python3 ${FMESHDRIVER} -n 20 20 --outDir ${CMAKE_CURRENT_BINARY_DIR}/full -s ${STENCILVAL} --bounds -5.0 5.0 -5.0 5.0 --periodic false")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Mesh generation failed")
else()
  message("Mesh generation succeeded!")
endif()

#------------------
# make sample mesh
#------------------
set(CMD "python3 ${SMESHDRIVER} --outdir ${CMAKE_CURRENT_BINARY_DIR} --fullMeshDir ${CMAKE_CURRENT_BINARY_DIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Mesh generation failed")
else()
  message("Mesh generation succeeded!")
endif()


#------------------
# run exe
#------------------
execute_process(COMMAND ${EXENAME} RESULT_VARIABLE CMD_RESULT)
if(RES)
  message(FATAL_ERROR "run failed")
else()
  message("run succeeded!")
endif()