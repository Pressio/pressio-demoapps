
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../../meshing_scripts)

add_subdirectory(mesh_s3_5x5x5_periodic)
add_subdirectory(mesh_s3_5x5x5_nonperiodic)

add_subdirectory(mesh_s5_5x5x5_periodic)
add_subdirectory(mesh_s5_5x5x5_nonperiodic)