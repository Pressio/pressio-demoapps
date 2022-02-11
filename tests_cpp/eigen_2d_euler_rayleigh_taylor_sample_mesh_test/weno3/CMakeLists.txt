
set(testname eigen_2d_euler_rayleigh_taylor_sample_mesh_weno3_test)
set(exename  ${testname}_exe)

configure_file(../../sample_mesh_compare.py sample_mesh_compare.py COPYONLY)
configure_file(../sample_mesh_gids.dat sample_mesh_gids.dat COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO3)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DFMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DSMESHDRIVER=${MESHSRC}/create_sample_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=5
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)
