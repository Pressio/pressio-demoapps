
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

add_subdirectory(mesh_generation_full_1d)

add_subdirectory(mesh_generation_full_2d)
add_subdirectory(mesh_generation_sample_2d)

add_subdirectory(mesh_generation_full_3d)
