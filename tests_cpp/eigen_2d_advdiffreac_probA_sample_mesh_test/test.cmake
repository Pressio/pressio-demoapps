include(FindUnixCommands)

#------------------
# compute full mesh
#------------------
set(CMD "python3 ${FMESHDRIVER} -n 13 13 --outDir ${OUTDIR} -s ${STENCILVAL} --bounds 0.0 1.0 0.0 1.0")
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

#------------------
# move files
#------------------
set(CMD "mkdir ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/state.txt ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/velo.txt ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/jacobian.txt ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/connectivity.dat ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/coordinates.dat ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

set(CMD "mv ${OUTDIR}/info.dat ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)

#------------------
# generate sample mesh
#------------------
set(CMD "python3 ${SMESHDRIVER} --outdir ${OUTDIR} --fullMeshDir ${OUTDIR}/full")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Mesh generation failed")
else()
  message("Mesh generation succeeded!")
endif()

#------------------
# run exe on sample mesh
#------------------
execute_process(COMMAND ${EXENAME} RESULT_VARIABLE CMD_RESULT)
if(RES)
  message(FATAL_ERROR "run failed")
else()
  message("run succeeded!")
endif()

#------------------
# check things are right
#------------------
set(CMD "python3 sample_mesh_compare.py --numdofspercell 1")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "comparison failed")
else()
  message("comparison succeeded!")
endif()
