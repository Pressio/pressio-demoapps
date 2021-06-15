
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../../meshing_scripts)

add_subdirectory(mesh_s3_5x1_periodic)
add_subdirectory(mesh_s3_5x5_periodic)
add_subdirectory(mesh_s3_5x6_periodic)
add_subdirectory(mesh_s3_6x5_periodic)
add_subdirectory(mesh_s5_8x8_periodic)
add_subdirectory(mesh_s5_7x8_periodic)
add_subdirectory(mesh_s5_8x7_periodic)
add_subdirectory(mesh_s7_8x8_periodic)
add_subdirectory(mesh_s7_7x8_periodic)
add_subdirectory(mesh_s7_8x7_periodic)

add_subdirectory(mesh_s3_5x6_nonperiodic)
