
set(testname eigen_1d_linear_adv_default_velo_sample_mesh_testA_weno5)
set(exename  ${testname}_exe)

configure_file(../../sample_mesh_compare.py sample_mesh_compare.py COPYONLY)
configure_file(../sample_mesh_gids.dat sample_mesh_gids.dat COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO5)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DFMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DSMESHDRIVER=${MESHSRC}/create_sample_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)
