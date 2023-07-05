include(FindUnixCommands)

# note that 24,32 is on purpose different to check that in 2d we do things right
set(CMD "python3 ${MESHDRIVER} -n 17 19 --outDir ${OUTDIR} -s ${STENCILVAL} --bounds -1.25 1.25 -1.25 1.25 --periodic x y")
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
