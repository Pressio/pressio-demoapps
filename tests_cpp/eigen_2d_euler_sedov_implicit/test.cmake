include(FindUnixCommands)

set(CMD "python3 ${MESHDRIVER} -n 18 18 --outDir ${OUTDIR} -s ${STENCILVAL} --bounds -0.5 0.5 -0.5 0.5")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Mesh generation failed")
else()
  message("Mesh generation succeeded!")
endif()

execute_process(COMMAND ${EXENAME} RESULT_VARIABLE CMD_RESULT)
if(RES)
  message(FATAL_ERROR "run failed")
else()
  message("run succeeded!")
endif()

set(CMD "python3 compare.py")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "comparison failed")
else()
  message("comparison succeeded!")
endif()
