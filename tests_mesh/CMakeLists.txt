
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../tpls/eigen3
  ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_subdirectory(mesh_generation_full_1d)
add_subdirectory(mesh_generation_full_2d)
add_subdirectory(mesh_generation_full_3d)
add_subdirectory(mesh_generation_sample_2d)
add_subdirectory(mesh_boundary_cells_tagging_2d)
