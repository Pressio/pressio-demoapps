include(FindUnixCommands)

set(MYLIST "10;20;40;80")
foreach(nx ${MYLIST})

  set(CMD "python3 ${MESHDRIVER} -n ${nx} 1 --outDir ${OUTDIR} -s ${SSTEN} --bounds -1.0 1.0 --periodic true")
  execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
  if(RES)
    message(FATAL_ERROR "Mesh generation failed")
  else()
    message("Mesh generation succeeded!")
  endif()

  execute_process(COMMAND ${EXENAME} RESULT_VARIABLE CMD_RESULT)
  if(RES)
    message(FATAL_ERROR "run with ${nx} failed")
  else()
    message("run with ${nx} succeeded!")
  endif()

  set(CMD "mv ${TESTNAME}_solution.bin ${TESTNAME}_sol${nx}.bin")
  execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
  set(CMD "mv coordinates.dat coords${nx}.dat")
  execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
endforeach()

# check convergence
set(CMD "python3 conv.py")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "test failed")
else()
  message("test succeeded!")
endif()
