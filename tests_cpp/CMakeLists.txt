
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

add_subdirectory(ccuMeshObj)
add_subdirectory(ccuMeshObjSpecialCellsLists)
add_subdirectory(weno5)

### 1d problems ###
add_subdirectory(1d_linear_adv_convergence)
add_subdirectory(1d_euler_convergence)
add_subdirectory(1d_sod)
add_subdirectory(1d_lax)

### 2d problems ###
add_subdirectory(2d_smooth)
add_subdirectory(2d_riemann)
add_subdirectory(2d_sedov)
add_subdirectory(2d_sedov_symmetry)
add_subdirectory(2d_euler_sample_mesh_test1)

### 3d problems ###
add_subdirectory(3d_smooth)
add_subdirectory(3d_sedov_symmetry)
