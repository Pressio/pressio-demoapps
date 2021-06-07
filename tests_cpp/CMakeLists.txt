
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

add_subdirectory(meshing)
add_subdirectory(weno5)
add_subdirectory(linadv1d_convergence)
add_subdirectory(euler1d_convergence)
add_subdirectory(sod1d)

add_subdirectory(sedov2d)
add_subdirectory(riemann2d)
