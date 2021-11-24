
#add_custom_command(COMMAND

# set(CMD "python ${MESHSRC}/create_full_mesh.py -n 25 25 --outDir ${CMAKE_CURRENT_BINARY_DIR} -s 7 --bounds -1.0 1.0 -1.0 1.0 --periodic true")
# add_custom_command(OUTPUT toni.txt
#   COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/gigi.txt

#   COMMENT "Creating "
#   )

# execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
# if(RES)
#   message(FATAL_ERROR "Mesh generation failed")
# else()
#   message("Mesh generation succeeded!")
# endif()

set(exename  eigen_2d_euler_periodic_firstorder_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

set(exename  eigen_2d_euler_periodic_weno3_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO3)

set(exename  eigen_2d_euler_periodic_weno5_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO5)
