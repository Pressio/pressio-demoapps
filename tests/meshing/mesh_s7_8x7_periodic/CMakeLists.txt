
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../compare_mesh_files.py compare_mesh_files.py COPYONLY)

configure_file(connectivity_gold.dat connectivity_gold.dat COPYONLY)
configure_file(coordinates_gold.dat  coordinates_gold.dat  COPYONLY)
configure_file(info_gold.dat	     info_gold.dat	   COPYONLY)

add_test(NAME full_mesh_8x7_s7_periodic
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-Dnx=8
-Dny=7
-Dss=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test2d.cmake
)
