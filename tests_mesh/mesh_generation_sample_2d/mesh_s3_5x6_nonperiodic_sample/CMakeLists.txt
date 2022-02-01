
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../compare_mesh_files.py compare_mesh_files.py COPYONLY)

configure_file(connectivity.dat connectivity_gold.dat COPYONLY)
configure_file(coordinates.dat  coordinates_gold.dat  COPYONLY)
configure_file(info.dat	     	info_gold.dat	   COPYONLY)
configure_file(stencil_mesh_gids.dat stencil_mesh_gids_gold.dat	 COPYONLY)
configure_file(sample_mesh_gids.dat sample_mesh_gids.dat	 COPYONLY)

add_test(NAME full_mesh_5x6_s3_nonperiodic_sample
COMMAND ${CMAKE_COMMAND}
-DFMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DSMESHDRIVER=${MESHSRC}/create_sample_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-Dnx=5
-Dny=6
-Dss=3
-DPER=none
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)
