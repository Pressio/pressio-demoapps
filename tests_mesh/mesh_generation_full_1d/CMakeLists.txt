
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../../meshing_scripts)

add_subdirectory(mesh_s3_5x1_periodic)
add_subdirectory(mesh_s5_5x1_periodic)
add_subdirectory(mesh_s7_5x1_periodic)

add_subdirectory(mesh_s3_5x1_nonperiodic)
add_subdirectory(mesh_s5_5x1_nonperiodic)
add_subdirectory(mesh_s7_5x1_nonperiodic)
