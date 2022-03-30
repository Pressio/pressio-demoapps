include(FindUnixCommands)

# generate the mesh
if(${PER} EQUAL 0)
  set(CMD "python3 ${MESHDRIVER} -n ${nx} 1 --outDir ${OUTDIR} -s ${ss} --bounds 0.0 1.0")
else()
  set(CMD "python3 ${MESHDRIVER} -n ${nx} 1 --outDir ${OUTDIR} -s ${ss} --bounds 0.0 1.0 --periodic ${PER}")
endif()
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
