include(FindUnixCommands)

# generate the mesh
set(CMD "python3 ${MESHDRIVER} -n ${nx} ${ny} --outDir ${OUTDIR} -s ${ss} --bounds 0.0 1.0 0.0 1.0 --periodic ${PER}")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Mesh generation failed")
else()
  message("Mesh generation succeeded!")
endif()

# use python script to make comparison
set(CMD "python3 compare_mesh_files.py 1e-10")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff failed")
else()
  message("Mesh files diff succeeded!")
endif()
