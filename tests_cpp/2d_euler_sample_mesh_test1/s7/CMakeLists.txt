
add_executable(veloSmExeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

configure_file(../compare.py compare.py COPYONLY)
configure_file(sample_mesh_gids.dat sample_mesh_gids.dat COPYONLY)

add_test(NAME euler2d_s7_samplemesh_test1
COMMAND ${CMAKE_COMMAND}
-DFMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DSMESHDRIVER=${MESHSRC}/create_sample_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:veloSmExeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)
